import adsk.core, adsk.fusion, traceback
import math, time
from collections import namedtuple
import pickle

COMMAND_ID = "trochoid_decelerator"
ID_NECESSARY_TAB = "necessary_tab_"
ID_OPTIONAL_TAB  = "optional_tab_"
ID_OPT_CH_GROUP    = ID_OPTIONAL_TAB + "centor_hole_group_"
ID_OPT_TGTOD_GROUP = ID_OPTIONAL_TAB + "trochoidal_gear_to_output_disk_group_"
ID_TV = "realtime_view"
ID_NES_IMG   = ID_NECESSARY_TAB + "description_image"
ID_NES_RR    = ID_NECESSARY_TAB + "reducation_ratio_"
ID_NES_EA    = ID_NECESSARY_TAB + "eccentric_amount_"
ID_NES_RGPD  = ID_NECESSARY_TAB + "ring_gear_pin_diameter_"
ID_NES_RGPPD = ID_NECESSARY_TAB + "ring_gear_pin_pitch_diameter_"
ID_NES_CGPN  = ID_NECESSARY_TAB + "trochoidal_gear_plot_num_"
ID_NES_MPA   = ID_NECESSARY_TAB + "minimum_Pressure_angle"
ID_OPT_IMG    = ID_OPTIONAL_TAB + "description_image"
ID_OPT_CGH_DR = ID_OPT_CH_GROUP + "draw_centor_hole"
ID_OPT_CGH_D  = ID_OPT_CH_GROUP + "centor_hole_diameter"
ID_OPT_DR_CAH  = ID_OPT_TGTOD_GROUP + "draw_trochoid_around_hole"
ID_OPT_DR_DP   = ID_OPT_TGTOD_GROUP + "draw_disk_pin"
ID_OPT_CHOTGOD = ID_OPT_TGTOD_GROUP + "choise_trochoidal_gear_or_output_disk"
ID_OPT_CHOTGOD_AN  = ID_OPT_TGTOD_GROUP + "around_hole_num"
ID_OPT_CHOTGOD_AD  = ID_OPT_TGTOD_GROUP + "around_hole_diameter"
ID_OPT_CHOTGOD_APD = ID_OPT_TGTOD_GROUP + "around_hole_position_diameter"
ID_OPT_CHOTGOD_ON  = ID_OPT_TGTOD_GROUP + "output_disk_num"
ID_OPT_CHOTGOD_OD  = ID_OPT_TGTOD_GROUP + "output_disk_diameter"
ID_OPT_CHOTGOD_OPD = ID_OPT_TGTOD_GROUP + "output_disk_position_diameter"

USER_CHANGEABLE_ID = [  ID_NES_RR, ID_NES_EA, ID_NES_RGPD, ID_NES_RGPPD, ID_NES_CGPN,
                        ID_OPT_CGH_DR, ID_OPT_CGH_D,
                        ID_OPT_DR_CAH, ID_OPT_CHOTGOD_AN, ID_OPT_CHOTGOD_AD, ID_OPT_CHOTGOD_APD,
                        ID_OPT_DR_DP, ID_OPT_CHOTGOD_ON, ID_OPT_CHOTGOD_OD, ID_OPT_CHOTGOD_OPD
                     ]

COMMAND_NAME = 'create cyclo reducer'
COMMAND_DESCRIPTION = 'サイクロ減速機用曲線作成スクリプト'

def compositeSimpson(func,upper,lower, splitNum):
    splitNum=int(splitNum)
    if splitNum&0b1:
        splitNum += 1
    h = (upper-lower)/splitNum
    ysum = func(lower) + 4*func(lower+h) + func(upper)
    for i in range(2,splitNum)[::2]:
        ysum += 2*func(lower+i*h) + 4*func(lower+(i+1)*h)
    return h/3*(ysum)

def bisectionMethod(func, upper, lower, maxError, maxCalcTimes=100):
    maxError=abs(maxError)
    calcTimes=0
    while True:
        calcTimes+=1
        x = (upper+lower)/2.0
        if (0.0 < func(x)*func(upper)):
            upper=x
        else:
            lower=x
        if (upper-lower<=maxError):
            return x
        elif (calcTimes==maxCalcTimes):
            return x

def numericalAnalysis(func, initialValue, maxError, maxCalcTimes=100):
    maxError = abs(maxError)
    calcTimes=0
    x = initialValue
    dfx = lambda x: (func(x*1.000001)-func(x*0.999999)) / (x*0.000002)
    while True:
        calcTimes+=1
        xn = x-func(x)/dfx(x)
        if abs(xn-x)<=maxError:
            return xn
        elif calcTimes==maxCalcTimes:
            return xn
        x = xn

class CycloidalReducer():
    def __init__(self, ringPinNum, ringPinRadius, ringPinPitchRadius, eccentricAmount):
        if(ringPinNum<2) or (ringPinRadius<=0) or (ringPinPitchRadius<=0) or (eccentricAmount<=0):
            raise ValueError("invalid argument")

        self.ringPinNum              = ringPinNum
        self.ringPinRadius           = ringPinRadius
        self.ringPinPitchRadius      = ringPinPitchRadius
        self.trochoidalGearThoothNum = ringPinNum-1
        self.eccentricAmount         = eccentricAmount
        self.reducationRatio = self.trochoidalGearThoothNum / (self.ringPinNum - self.trochoidalGearThoothNum)
        self.rm = self.ringPinPitchRadius/(self.reducationRatio+1)
        self.rc = self.rm*self.reducationRatio
        self.rd = self.eccentricAmount
        self.d  = self.ringPinRadius

    def hasSingularPoint(self):
        if(self.eccentricAmount*self.ringPinNum >= self.ringPinPitchRadius):
            return True
        if self.getMaxOffset() < self.ringPinRadius:
            return True
        return False

    def fxa(self, p):
        return (self.rc+self.rm)*math.cos(p) - self.rd*math.cos((self.rc+self.rm)/self.rm*p)
    def fya(self, p):
        return (self.rc+self.rm)*math.sin(p) - self.rd*math.sin((self.rc+self.rm)/self.rm*p)
    def dfxa(self, p):
        return -(self.rc+self.rm)*math.sin(p) + ((self.rc+self.rm)/self.rm)*self.rd*math.sin((self.rc+self.rm)/self.rm*p)
    def dfya(self, p):
        return (self.rc+self.rm)*math.cos(p) - ((self.rc+self.rm)/self.rm)*self.rd*math.cos((self.rc+self.rm)/self.rm*p)
    def ddfxa(self, p):
        return -(self.rc+self.rm)*math.cos(p) + ((self.rc+self.rm)/self.rm)**2 * self.rd*math.cos((self.rc+self.rm)/self.rm*p)
    def ddfya(self, p):
        return -(self.rc+self.rm)*math.sin(p) + ((self.rc+self.rm)/self.rm)**2 * self.rd*math.sin((self.rc+self.rm)/self.rm*p)
    def dddfxa(self, p):
        return +(self.rc+self.rm)*math.sin(p) - ((self.rc+self.rm)/self.rm)**3 * self.rd*math.sin((self.rc+self.rm)/self.rm*p)
    def dddfya(self, p):
        return -(self.rc+self.rm)*math.cos(p) + ((self.rc+self.rm)/self.rm)**3 * self.rd*math.cos((self.rc+self.rm)/self.rm*p)

    def fxp(self, p):
        dxa, dya = self.dfxa(p), self.dfya(p)
        return self.fxa(p) - self.d*dya / math.sqrt(dxa**2+dya**2)
    def fyp(self, p):
        dxa, dya = self.dfxa(p), self.dfya(p)
        return self.fya(p) + self.d*dxa / math.sqrt(dxa**2+dya**2)
    def dfxp(self, p):
        dxa,  dya  = self.dfxa(p),  self.dfya(p)
        ddxa, ddya = self.ddfxa(p), self.ddfya(p)
        D = self.d
        return dxa * (1 + D*(-dxa*ddya + dya*ddxa)/(dxa**2+dya**2)**(3/2) )
    def dfyp(self, p):
        dxa,  dya  = self.dfxa(p),  self.dfya(p)
        ddxa, ddya = self.ddfxa(p), self.ddfya(p)
        D    = self.d
        return dya * (1 + D*(dya*ddxa - dxa*ddya)/(dxa**2+dya**2)**(3/2) )
    def ddfxp(self, p):
        dxa,   dya   = self.dfxa(p),   self.dfya(p)
        ddxa,  ddya  = self.ddfxa(p),  self.ddfya(p)
        dddxa, dddya = self.dddfxa(p), self.dddfya(p)
        D = self.d
        w = dxa**2+dya**2
        return ddxa   * 1 -dddya * D*w**-0.5 -ddya  * -2*D*(dxa*ddxa+dya*ddya)*w**-1.5 -dya   * -D*((ddxa**2 + dxa*dddxa + ddya**2 + dya*dddya)*w**-1.5 -3*(dxa*ddxa+dya*ddya)**2*w**-2.5)
    def ddfyp(self, p):
        dxa,   dya   = self.dfxa(p),   self.dfya(p)
        ddxa,  ddya  = self.ddfxa(p),  self.ddfya(p)
        dddxa, dddya = self.dddfxa(p), self.dddfya(p)
        D = self.d
        w = dxa**2+dya**2
        return +ddya  * 1 +dddxa * D*w**-0.5 +ddxa  * -2*D*(dxa*ddxa+dya*ddya)*w**-1.5 +dxa   * -D*((ddxa**2 + dxa*dddxa + ddya**2 + dya*dddya)*w**-1.5 -3*(dxa*ddxa+dya*ddya)**2*w**-2.5)

    def fa(self, p):
        a = math.atan2(self.dfyp(p),self.dfxp(p)) - math.atan2(self.fyp(p),self.fxp(p))
        return a
    def dfa(self, p):
        xp, yp  = self.fxp(p), self.fyp(p)
        dxp, dyp  = self.dfxp(p), self.dfyp(p)
        ddxp, ddyp = self.ddfxp(p), self.ddfyp(p)
        return (dxp*ddyp - ddxp*dyp)/(dxp**2+dyp**2) - (xp*dyp - dxp*yp)/(xp**2+yp**2)
    def getMinimumPressureAngle(self):
        lastPOneThooth = 2*math.pi/self.trochoidalGearThoothNum
        maxError=0.00001
        minP = bisectionMethod(self.dfa, lastPOneThooth/2, 0, maxError)
        minAngle = self.fa(minP)
        return minAngle if minAngle<=math.pi/2.0 else math.pi-minAngle

    def fcr(self, p):
        dxa, dya  = self.dfxa(p), self.dfya(p)
        ddxa, ddya = self.ddfxa(p), self.ddfya(p)
        numerator = (dxa**2 + dya**2)**1.5
        denominator = dxa*ddya - dya*ddxa
        return numerator / denominator

    def getMaxOffset(self):
        (rc, rm, rd) = (self.rc, self.rm, self.rd)
        inacos = (2*rc*rd**2 - rc*rm**2 + rd**2*rm + rm**3)/(rd*rm*(rc+2*rm))
        if abs(inacos) <= 1:
            s = [self.fcr(rm/rc*math.pi), self.fcr(rm/rc*math.acos(inacos))]
            return min(s)
        else:
            return self.fcr(rm/rc*math.pi)

    def getPerimeter(self, upper, lower, splitNum=1000):
        dfl = lambda p: math.sqrt(self.dfxp(p)**2 + self.dfyp(p)**2)
        return compositeSimpson(dfl, upper, lower, splitNum)
    def getConstDistancePoint(self, currentP, distance, upperP):
        f = lambda p: self.getPerimeter(p, currentP, 100)-distance
        return bisectionMethod(f, upperP, currentP, self.pointError)

    def getTrochoidPoints(self, pointNum, shift=False):
        centor = [self.eccentricAmount,0] if shift else [0,0]
        lastPOneThooth = 2*math.pi/self.trochoidalGearThoothNum
        pOneThooth = [i*lastPOneThooth/pointNum for i in range(pointNum)]

        pAllThooth = []
        for i in range(self.trochoidalGearThoothNum):
            pAllThooth += [i*lastPOneThooth+p for p in pOneThooth]

        points=[]
        for p in pAllThooth:
            points.append([self.fxa(p)+centor[0], self.fya(p)+centor[1]])
        return (points, centor)

    def getTrochoidParallelCurvePoints(self, pointNum, shift=True):
        centor = [self.eccentricAmount,0] if shift else [0,0]
        lastPOneThooth = 2*math.pi/self.trochoidalGearThoothNum

        perimeterOneThooth = self.getPerimeter(lastPOneThooth, 0, 1000)/pointNum
        self.pointError = perimeterOneThooth/pointNum/1000000

        pOneThooth=[0]
        for i in range(pointNum)[1:]:
            px = pOneThooth[-1]
            pOneThooth.append( self.getConstDistancePoint(px, perimeterOneThooth, lastPOneThooth) )

        pAllThooth = []
        for i in range(self.trochoidalGearThoothNum):
            pAllThooth += [i*lastPOneThooth+p for p in pOneThooth]

        points=[]
        for p in pAllThooth:
            points.append([self.fxp(p)+centor[0], self.fyp(p)+centor[1]])
        return (points, centor)

    def getOutpinPoints(self):
        points = []
        for i in range(self.ringPinNum):
            theta = 2*math.pi * (i/self.ringPinNum)
            x = self.ringPinPitchRadius*math.cos(theta)
            y = self.ringPinPitchRadius*math.sin(theta)
            points.append([x,y])
        return (points, self.ringPinRadius)

class DrawCycloReducer():
    def __init__(self, inputs):
        drawingParam = inputsToParameter(inputs)

        design = _app.activeProduct
        activeComp = design.activeOccurrence.component if design.activeOccurrence else  design.rootComponent
        occTrochoidalGear = activeComp.occurrences.addNewComponent(adsk.core.Matrix3D.create())
        compReducer = occTrochoidalGear.component
        compReducer.name = "Cycloidal reducer"

        self.cycoroidDecelerator = CycloidalReducer(int(drawingParam.ringPinNum), drawingParam.ringPinDia/2.0,
                                                    drawingParam.ringPinPitchDia/2.0,  drawingParam.eccentricAmount)

        trochoidSketch = compReducer.sketches.add(compReducer.xYConstructionPlane)
        ringPinSketch = compReducer.sketches.add(compReducer.xYConstructionPlane)
        skts = [trochoidSketch, ringPinSketch]

        rr   = int(drawingParam.ringPinNum-1)
        ea   = _unitsMgr.convert(drawingParam.eccentricAmount, _unitsMgr.internalUnits, _unitsMgr.defaultLengthUnits)
        rpd  = _unitsMgr.convert(drawingParam.ringPinDia,      _unitsMgr.internalUnits, _unitsMgr.defaultLengthUnits)
        rppd = _unitsMgr.convert(drawingParam.ringPinPitchDia, _unitsMgr.internalUnits, _unitsMgr.defaultLengthUnits)
        eaString   = "{:.3g}".format(ea)   + _unitsMgr.defaultLengthUnits
        rpdString  = "{:.3g}".format(rpd)  + _unitsMgr.defaultLengthUnits
        rppdString = "{:.3g}".format(rppd) + _unitsMgr.defaultLengthUnits
        trochoidSketch.name = "Cycloidal gear"+"(rr:"+str(rr)+" ea:"+eaString+")"
        ringPinSketch.name = "Ring pins"+"(rpd:"+rpdString+" rppd:"+rppdString+")"

        if drawingParam.isDrawOutputDiskPin:
            outputDiskSketch = compReducer.sketches.add(compReducer.xYConstructionPlane)
            outputDiskSketch.name = "Output disk pin"
            skts.append(outputDiskSketch)

        for skt in skts:
            skt.isComputeDeferred = True

        self.createTrochoidalGear(trochoidSketch, drawingParam)
        self.createRingGear(ringPinSketch, drawingParam)
        if drawingParam.isDrawCentorHole:
            self.createTrochoidalGearCentorHole(trochoidSketch, drawingParam)
        if drawingParam.isDrawAroundHole:
            self.createTrochoidalGearAroundHole(trochoidSketch, drawingParam)
        if drawingParam.isDrawOutputDiskPin:
            self.createOutputDisk(outputDiskSketch, drawingParam)

        for skt in skts:
            skt.isComputeDeferred = False

    def createTrochoidalGear(self, sketch, drawingParam):
        sketchOriginPoint = sketch.originPoint
        z=0
        (trochoidParallelPoints, trochoidalGearCentor) = self.cycoroidDecelerator.getTrochoidParallelCurvePoints(drawingParam.plotDotNum)

        self.trochoidCentorPoint2D = sketch.sketchPoints.add( adsk.core.Point3D.create(trochoidalGearCentor[0], trochoidalGearCentor[1], z) )
        self.distanceDimentionEasy(sketch, self.trochoidCentorPoint2D, sketchOriginPoint)

        splinePoints = adsk.core.ObjectCollection.create()
        for (x,y) in trochoidParallelPoints:
            splinePoints.add(adsk.core.Point3D.create(x, y, z))
        trochoidCurve = sketch.sketchCurves.sketchFittedSplines.add(splinePoints)
        trochoidCurve.isClosed = True
        trochoidCurve.isFixed  = True

        if False:
            (trochoidPoints, trochoidalGearCentor) = self.cycoroidDecelerator.getTrochoidPoints(drawingParam.plotDotNum, True)
            splinePoints2 = adsk.core.ObjectCollection.create()
            for (xa,ya) in trochoidPoints:
                splinePoints2.add(adsk.core.Point3D.create(xa, ya, z))
            trochoidCurve2 = sketch.sketchCurves.sketchFittedSplines.add(splinePoints2)
            trochoidCurve2.isClosed = True
            trochoidCurve2.isFixed  = True
            trochoidCurve2.isConstruction = True

    def createTrochoidalGearAroundHole(self, sketch, drawingParam):
        sketchOriginPoint = sketch.originPoint
        n = drawingParam.troGearAroundHoleNum
        r = drawingParam.troGearAroundHoleDia/2.0
        pd = drawingParam.troGearAroundHolePosDia
        pdn = drawingParam.plotDotNum
        z=0

        (trochoidalGearPoints, trochoidalGearCentor) = self.cycoroidDecelerator.getTrochoidPoints(pdn, True)

        trochoidCentorPoint2D = sketch.sketchPoints.add( adsk.core.Point3D.create(trochoidalGearCentor[0], trochoidalGearCentor[1], z) )
        self.distanceDimentionEasy(sketch, trochoidCentorPoint2D, sketchOriginPoint)
        fx = lambda theta : pd/2.0 * math.cos(theta) + trochoidalGearCentor[0]
        fy = lambda theta : pd/2.0 * math.sin(theta) + trochoidalGearCentor[1]
        z=0
        firstCircle = sketch.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(fx(0), fy(0), z), r)
        self.distanceDimentionEasy(sketch, trochoidCentorPoint2D, firstCircle.centerSketchPoint)
        self.diameterDimentionEasy(sketch, firstCircle)
        firstLine = sketch.sketchCurves.sketchLines.addByTwoPoints(trochoidCentorPoint2D, firstCircle.centerSketchPoint)
        firstLine.isConstruction=True

        beforeLine = firstLine
        beforeAngleDim = None
        for i in range(n)[1:]:
            theta = (float(i)/n)*2*math.pi
            circle = sketch.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(fx(theta), fy(theta), z), r)
            sketch.geometricConstraints.addEqual(firstCircle, circle)
            line = sketch.sketchCurves.sketchLines.addByTwoPoints(trochoidCentorPoint2D, circle.centerSketchPoint)
            line.isConstruction=True
            sketch.geometricConstraints.addEqual(firstLine, line)

            angleDim = self.angleDimentionEasy(sketch, beforeLine, line)
            if beforeAngleDim:
                angleDim.parameter.expression = beforeAngleDim.parameter.name

            beforeLine = line
            beforeAngleDim = angleDim

    def createTrochoidalGearCentorHole(self, sketch, drawingParam):
        sketchOriginPoint = sketch.originPoint
        hd = drawingParam.troGearCentorHoleDia
        pdn = drawingParam.plotDotNum
        z=0

        (trochoidalGearPoints, trochoidalGearCentor) = self.cycoroidDecelerator.getTrochoidPoints(pdn, True)

        p = adsk.core.Point3D.create(trochoidalGearCentor[0], trochoidalGearCentor[1], z)
        circle = sketch.sketchCurves.sketchCircles.addByCenterRadius(p, hd/2.0)
        self.distanceDimentionEasy(sketch, sketchOriginPoint, circle.centerSketchPoint)
        return self.diameterDimentionEasy(sketch, circle)

    def createRingGear(self, sketch, drawingParam):
        sketchOriginPoint = sketch.originPoint
        comp = sketch.parentComponent
        z=0
        (points, radius) = self.cycoroidDecelerator.getOutpinPoints()

        firstThoothXY = points[0]
        firstCircle = sketch.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(firstThoothXY[0], firstThoothXY[1], z), radius)
        self.distanceDimentionEasy(sketch, sketchOriginPoint, firstCircle.centerSketchPoint)
        self.diameterDimentionEasy(sketch, firstCircle)
        firstLine = sketch.sketchCurves.sketchLines.addByTwoPoints(sketchOriginPoint, firstCircle.centerSketchPoint)
        firstLine.isConstruction=True

        beforeLine = firstLine
        beforeAngleDim = None
        for p in points[1:]:
            circle = sketch.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(p[0], p[1], z), radius)
            sketch.geometricConstraints.addEqual(firstCircle, circle)

            line = sketch.sketchCurves.sketchLines.addByTwoPoints(sketchOriginPoint, circle.centerSketchPoint)
            line.isConstruction=True
            sketch.geometricConstraints.addEqual(firstLine, line)

            angleDim = self.angleDimentionEasy(sketch, beforeLine, line)
            if beforeAngleDim:
                angleDim.parameter.expression = beforeAngleDim.parameter.name

            beforeLine = line
            beforeAngleDim = angleDim

    def createOutputDisk(self, sketch, drawingParam):
        sketchOriginPoint = sketch.originPoint
        z=0
        (points, radius) = self.cycoroidDecelerator.getOutpinPoints()

        n              = drawingParam.outDiskPinNum
        positionRadius = drawingParam.outDiskPinPosDia/2.0
        holeRadius     = drawingParam.outDiskPinDia/2.0
        fx = lambda theta : positionRadius * math.cos(theta)
        fy = lambda theta : positionRadius * math.sin(theta)

        firstCircle = sketch.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(fx(0), fy(0), z), holeRadius)
        self.distanceDimentionEasy(sketch, sketchOriginPoint, firstCircle.centerSketchPoint)
        self.diameterDimentionEasy(sketch, firstCircle)
        firstLine = sketch.sketchCurves.sketchLines.addByTwoPoints(sketchOriginPoint, firstCircle.centerSketchPoint)
        firstLine.isConstruction=True

        beforeLine = firstLine
        beforeAngleDim = None
        for i in range(n)[1:]:
            theta = i/n*2*math.pi
            circle = sketch.sketchCurves.sketchCircles.addByCenterRadius(adsk.core.Point3D.create(fx(theta), fy(theta), z), holeRadius)
            sketch.geometricConstraints.addEqual(firstCircle, circle)
            line = sketch.sketchCurves.sketchLines.addByTwoPoints(sketchOriginPoint, circle.centerSketchPoint)
            line.isConstruction=True
            sketch.geometricConstraints.addEqual(firstLine, line)

            angleDim = self.angleDimentionEasy(sketch, beforeLine, line)
            if beforeAngleDim:
                angleDim.parameter.expression = beforeAngleDim.parameter.name

            beforeLine = line
            beforeAngleDim = angleDim

    def distanceDimentionEasy(self, sketch, sketchPoint1, sketchPoint2):
        z=0
        textPoint3D = adsk.core.Point3D.create((sketchPoint1.geometry.x+sketchPoint2.geometry.x)/2,
                                               (sketchPoint1.geometry.y+sketchPoint2.geometry.y)/2, z)
        sketch.sketchDimensions.addDistanceDimension(sketchPoint1, sketchPoint2,
                                                     adsk.fusion.DimensionOrientations.HorizontalDimensionOrientation,
                                                     textPoint3D)
        sketch.sketchDimensions.addDistanceDimension(sketchPoint1,sketchPoint2,
                                                     adsk.fusion.DimensionOrientations.VerticalDimensionOrientation,
                                                     textPoint3D)
    def diameterDimentionEasy(self, sketch, circle):
        z=0
        r = circle.radius
        textPoint3D = adsk.core.Point3D.create(circle.centerSketchPoint.geometry.x+r/2,
                                               circle.centerSketchPoint.geometry.y+r/2, z)
        sketch.sketchDimensions.addDiameterDimension(circle, textPoint3D)
    def angleDimentionEasy(self, sketch, line1, line2):
        z=0
        textPoint3D = adsk.core.Point3D.create((line1.geometry.startPoint.x+line1.geometry.endPoint.x+line2.geometry.startPoint.x+line2.geometry.endPoint.x)/3/2,
                                               (line1.geometry.startPoint.y+line1.geometry.endPoint.y+line2.geometry.startPoint.y+line2.geometry.endPoint.y)/3/2,
                                               (line1.geometry.startPoint.z+line1.geometry.endPoint.z+line2.geometry.startPoint.z+line2.geometry.endPoint.z)/3/2)
        return sketch.sketchDimensions.addAngularDimension(line1, line2, textPoint3D)

def inputsToParameter(commandInputs):
    drawingParam = namedtuple("DrawingParam",
                            ("ringPinNum", "ringPinDia", "ringPinPitchDia",
                                "eccentricAmount", "plotDotNum",
                                "troGearAroundHoleNum", "troGearAroundHoleDia", "troGearAroundHolePosDia",
                                "troGearCentorHoleDia",
                                "outDiskPinNum", "outDiskPinDia","outDiskPinPosDia"
                                "isDrawTrochoidalGear", "isDrawRingPin","isDrawCentorHole", "isDrawAroundHole","isDrawOutputDiskPin"
                                ))

    drawingParam.isDrawTrochoidalGear = True
    drawingParam.isDrawRingPin       = True
    drawingParam.isDrawCentorHole    = commandInputs.itemById(ID_OPT_CGH_DR).value
    drawingParam.isDrawAroundHole    = commandInputs.itemById(ID_OPT_DR_CAH).value
    drawingParam.isDrawOutputDiskPin = commandInputs.itemById(ID_OPT_DR_DP).value

    reducationRatioInput   = commandInputs.itemById(ID_NES_RR)
    eccentricAmountInput   = commandInputs.itemById(ID_NES_EA)
    ringPinDiaInput    = commandInputs.itemById(ID_NES_RGPD)
    ringPinPitchDiaInput = commandInputs.itemById(ID_NES_RGPPD)
    plotNumInput           = commandInputs.itemById(ID_NES_CGPN)

    drawingParam.ringPinNum    = int(reducationRatioInput.value)+1
    drawingParam.ringPinDia    = _unitsMgr.evaluateExpression(ringPinDiaInput.expression)
    drawingParam.ringPinPitchDia = _unitsMgr.evaluateExpression(ringPinPitchDiaInput.expression)
    drawingParam.eccentricAmount   = _unitsMgr.evaluateExpression(eccentricAmountInput.expression)
    drawingParam.plotDotNum        = int(plotNumInput.value)

    if drawingParam.isDrawCentorHole:
        troGearCentorHoleDiaInput    = commandInputs.itemById(ID_OPT_CGH_D)
        drawingParam.troGearCentorHoleDia = _unitsMgr.evaluateExpression(troGearCentorHoleDiaInput.expression)

    if drawingParam.isDrawAroundHole:
        troGearAroundHoleNumInput    = commandInputs.itemById(ID_OPT_CHOTGOD_AN)
        troGearAroundHoleDiaInput    = commandInputs.itemById(ID_OPT_CHOTGOD_AD)
        troGearAroundHolePosDiaInput = commandInputs.itemById(ID_OPT_CHOTGOD_APD)

        drawingParam.troGearAroundHoleNum    = int(troGearAroundHoleNumInput.value)
        drawingParam.troGearAroundHoleDia    = _unitsMgr.evaluateExpression(troGearAroundHoleDiaInput.expression)
        drawingParam.troGearAroundHolePosDia = _unitsMgr.evaluateExpression(troGearAroundHolePosDiaInput.expression)

    if drawingParam.isDrawOutputDiskPin:
        outDiskPinNumInput    = commandInputs.itemById(ID_OPT_CHOTGOD_ON)
        outDiskPinDiaInput    = commandInputs.itemById(ID_OPT_CHOTGOD_OD)
        outDiskPinPosDiaInput = commandInputs.itemById(ID_OPT_CHOTGOD_OPD)

        drawingParam.outDiskPinNum    = int(outDiskPinNumInput.value)
        drawingParam.outDiskPinDia    = _unitsMgr.evaluateExpression(outDiskPinDiaInput.expression)
        drawingParam.outDiskPinPosDia = _unitsMgr.evaluateExpression(outDiskPinPosDiaInput.expression)
    return drawingParam

def settingComandInputsItem(inputs):
    testViewInputs = inputs.addBoolValueInput(ID_TV, "Test view", False, "", False)
    testViewInputs.isFullWidth = True
    necessaryTabInput = inputs.addTabCommandInput(ID_NECESSARY_TAB, "Necessary param")
    necessaryTabChildInputs = necessaryTabInput.children
    necImageInputs = necessaryTabChildInputs.addImageCommandInput(ID_NES_IMG, "", "image/cyclo_nec.png")
    necImageInputs.isFullWidth = True
    reductionRatioInput = necessaryTabChildInputs.addIntegerSpinnerCommandInput(ID_NES_RR, 'Reduction ratio', 2, 99999, 1, 51)
    reductionRatioInput.tooltip = "ReductionRatio = RingPinNum-1 = cycloidalGear's tooth num"
    necessaryTabChildInputs.addValueInput(ID_NES_EA,   "Eccentric amount",        "mm", adsk.core.ValueInput.createByReal(1.5))
    necessaryTabChildInputs.addValueInput(ID_NES_RGPD, 'Ring pin diameter',       'mm', adsk.core.ValueInput.createByReal(6.0))
    necessaryTabChildInputs.addValueInput(ID_NES_RGPPD,'Ring pin pitch diameter', 'mm', adsk.core.ValueInput.createByReal(160.0))
    necessaryTabChildInputs.addIntegerSpinnerCommandInput(ID_NES_CGPN, "Cycloidal curve plot num per tooth", 2, 99999, 1, 6)
    necessaryTabChildInputs.addTextBoxCommandInput(ID_NES_MPA, "minimum Pressure angle", "-", 1, True)
    optionTabInput = inputs.addTabCommandInput(ID_OPTIONAL_TAB, "Optional param")
    optionTabChildInputs = optionTabInput.children
    optImageInputs = optionTabChildInputs.addImageCommandInput(ID_OPT_IMG, "", "image/cyclo_opt.png")
    optImageInputs.isFullWidth = True
    centorHoleGroup = optionTabChildInputs.addGroupCommandInput(ID_OPT_CH_GROUP, "Cycloidal gear center hole")
    centorHoleInputs = centorHoleGroup.children
    centorHoleInputs.addBoolValueInput(ID_OPT_CGH_DR, "Draw center hole", True, "", False)
    centorHoleInputs.addValueInput(ID_OPT_CGH_D, "Diameter", "mm", adsk.core.ValueInput.createByReal(55.0))
    trochoidToOutputGroup = optionTabChildInputs.addGroupCommandInput(ID_OPT_TGTOD_GROUP, "Cycloidal gear to output disk")
    trochoidToOutputInputs = trochoidToOutputGroup.children
    trochoidToOutputInputs.addBoolValueInput(ID_OPT_DR_CAH, "Draw around hole", True, "", False)
    trochoidToOutputInputs.addBoolValueInput(ID_OPT_DR_DP, "Draw output disk pin", True, "", False)
    holeOrPinSelectInputs = trochoidToOutputInputs.addDropDownCommandInput(ID_OPT_CHOTGOD, "Set about", adsk.core.DropDownStyles.LabeledIconDropDownStyle)
    holeOrPinSelectListItems = holeOrPinSelectInputs.listItems
    holeOrPinSelectListItems.add("Cycloidal gear hole", True)
    holeOrPinSelectListItems.add("Output disk pin", False)
    trochoidToOutputInputs.addIntegerSpinnerCommandInput(ID_OPT_CHOTGOD_AN, "Hole num", 2, 99999, 1, 10)
    trochoidToOutputInputs.addValueInput(ID_OPT_CHOTGOD_AD,  "Hole diameter",           "mm", adsk.core.ValueInput.createByReal(12.0))
    trochoidToOutputInputs.addValueInput(ID_OPT_CHOTGOD_APD, "Center to hole distance", "mm", adsk.core.ValueInput.createByReal(70.0))
    trochoidToOutputInputs.addIntegerSpinnerCommandInput(ID_OPT_CHOTGOD_ON, "Pin num", 2, 99999, 1, 10)
    trochoidToOutputInputs.addValueInput(ID_OPT_CHOTGOD_OD,  "Pin diameter",            "mm", adsk.core.ValueInput.createByReal(8.0))
    trochoidToOutputInputs.addValueInput(ID_OPT_CHOTGOD_OPD, "Center to pin distance",  "mm", adsk.core.ValueInput.createByReal(70.0))

handlers = []

class MyCommandCreatedHandler(adsk.core.CommandCreatedEventHandler):
    def __init__(self):
        super().__init__()
    def notify(self, args):
        try:
            cmd = args.command
            cmd.setDialogInitialSize(300,500)
            onExecute = MyCommandExecuteHandler()
            cmd.execute.add(onExecute)
            onDestroy = MyCommandDestroyHandler()
            cmd.destroy.add(onDestroy)
            onValidateInputs = MyCommandValidateInputsHandler()
            cmd.validateInputs.add(onValidateInputs)
            onExecutePreview = MyExecutePreviewHandler()
            cmd.executePreview.add(onExecutePreview)

            handlers.append(onExecute)
            handlers.append(onDestroy)
            handlers.append(onValidateInputs)
            handlers.append(onExecutePreview)

            inputs = cmd.commandInputs
            settingComandInputsItem(inputs)
            loadInputsValue(inputs, USER_CHANGEABLE_ID)

        except:
            if _ui:
                _ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

def saveInputsValues(commandInputs, ids, groupName="lastCommandInputValue"):
    design = _app.activeProduct
    attribs = design.attributes
    for id in ids:
        value = commandInputs.itemById(id).value
        binaryValue = pickle.dumps(value)
        strValue = binaryValue.hex()
        attribs.add(groupName, id, strValue)

def loadInputsValue(commandInputs, ids, groupName="lastCommandInputValue"):
    design = _app.activeProduct
    attribs = design.attributes
    if groupName not in attribs.groupNames:
        return False
    for id in ids:
        strValue = attribs.itemByName(groupName, id).value
        binaryValue = bytes.fromhex(strValue)
        value = pickle.loads(binaryValue)
        commandInputs.itemById(id).value = value
    return True

class MyCommandValidateInputsHandler(adsk.core.ValidateInputsEventHandler):
    def __init__(self):
        super().__init__()
    def notify(self, args):
        try:
            args.areInputsValid = True
            cmd = args.firingEvent.sender
            inputs = cmd.commandInputs
            param  = inputsToParameter(inputs)

            try:
                cr = CycloidalReducer(param.ringPinNum, param.ringPinDia/2.0, param.ringPinPitchDia/2.0, param.eccentricAmount)
                if cr.hasSingularPoint():
                    args.areInputsValid = False
            except ValueError as e:
                print(str(e))
                args.areInputsValid = False
            if (param.eccentricAmount <=0) or (param.ringPinPitchDia <=0) or (param.ringPinDia <=0) or (param.ringPinNum<2) or (param.eccentricAmount*param.ringPinNum >= (param.ringPinPitchDia/2.0)):
                args.areInputsValid = False

            if param.isDrawCentorHole is True:
                if param.troGearCentorHoleDia <=0:
                    args.areInputsValid = False
            if param.isDrawAroundHole is True:
                if param.troGearAroundHoleNum<=0 or param.troGearAroundHoleDia<=0 or param.troGearAroundHolePosDia<=0:
                    args.areInputsValid = False
            if param.isDrawOutputDiskPin is True:
                if param.outDiskPinNum<=0 or param.outDiskPinDia<=0 or param.outDiskPinPosDia<=0:
                    args.areInputsValid = False

            centorHoleDiaInput  = inputs.itemById(ID_OPT_CGH_D)
            centorHoleDiaInput.isEnabled = inputs.itemById(ID_OPT_CGH_DR).value

            trochoidOrOutputInput = inputs.itemById(ID_OPT_CHOTGOD)
            drawAroundHoleInput   = inputs.itemById(ID_OPT_DR_CAH)
            drawOutputPinInput    = inputs.itemById(ID_OPT_DR_DP)
            AroundHoleNumInput         = inputs.itemById(ID_OPT_CHOTGOD_AN)
            AroundHoleDiameterInput    = inputs.itemById(ID_OPT_CHOTGOD_AD)
            AroundHolePosDiameterInput = inputs.itemById(ID_OPT_CHOTGOD_APD)
            outputNumInput         = inputs.itemById(ID_OPT_CHOTGOD_ON)
            outputDiameterInput    = inputs.itemById(ID_OPT_CHOTGOD_OD)
            outputPosDiameterInput = inputs.itemById(ID_OPT_CHOTGOD_OPD)

            isDrawAroundHole = inputs.itemById(ID_OPT_DR_CAH).value
            isDrawOutputDisk = inputs.itemById(ID_OPT_DR_DP).value

            aroundInputList = [AroundHoleNumInput, AroundHoleDiameterInput, AroundHolePosDiameterInput]
            outputInputList = [outputNumInput, outputDiameterInput, outputPosDiameterInput]

            eccentricAmountInput = inputs.itemById(ID_NES_EA)
            eccentricAmountCm = _unitsMgr.evaluateExpression(eccentricAmountInput.expression)

            if (drawAroundHoleInput.value is False) and (drawOutputPinInput.value is False):
                trochoidOrOutputInput.isEnabled = False
                for i in aroundInputList + outputInputList:
                    i.isEnabled = False
            else:
                trochoidOrOutputInput.isEnabled = True
                for i in aroundInputList + outputInputList:
                    i.isEnabled = True

            aroundDiaCm  = _unitsMgr.evaluateExpression(AroundHoleDiameterInput.expression)
            outputDiaCm  = _unitsMgr.evaluateExpression(outputDiameterInput.expression)
            if trochoidOrOutputInput.selectedItem.index == 0:
                for i in aroundInputList:
                    i.isVisible = True
                for i in outputInputList:
                    i.isVisible = False
                outputNumInput.value              = AroundHoleNumInput.value
                outputDiameterInput.expression    = _unitsMgr.formatInternalValue(aroundDiaCm - 2*eccentricAmountCm)
                outputPosDiameterInput.expression = AroundHolePosDiameterInput.expression

            else:
                for i in aroundInputList:
                    i.isVisible = False
                for i in outputInputList:
                    i.isVisible = True
                AroundHoleNumInput.value              = outputNumInput.value
                AroundHoleDiameterInput.expression    = _unitsMgr.formatInternalValue(outputDiaCm + 2*eccentricAmountCm)
                AroundHolePosDiameterInput.expression = outputPosDiameterInput.expression

            if isDrawAroundHole:
                if (AroundHoleNumInput.value <= 0) or (aroundDiaCm <= 0) or (_unitsMgr.evaluateExpression(AroundHolePosDiameterInput.expression) <= 0):
                    args.areInputsValid = False
            if isDrawOutputDisk:
                if (outputNumInput.value <= 0) or (outputDiaCm <= 0) or (_unitsMgr.evaluateExpression(outputPosDiameterInput.expression) <= 0):
                    args.areInputsValid = False

        except:
            if _ui:
                _ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

class MyExecutePreviewHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()
    def notify(self, args):
        try:
            command = args.firingEvent.sender
            inputs = command.commandInputs

            PressureAngleInput = inputs.itemById(ID_NES_MPA)
            drawingParam = inputsToParameter(inputs)
            cycoroidDecelerator = CycloidalReducer(int(drawingParam.ringPinNum), drawingParam.ringPinDia/2.0,
                                                    drawingParam.ringPinPitchDia/2.0,  drawingParam.eccentricAmount)
            PressureAngle = cycoroidDecelerator.getMinimumPressureAngle()*180/math.pi
            PressureAngleInput.text = str( round(PressureAngle, 2))

            testView = inputs.itemById(ID_TV)
            if testView.value:
                testView.value = False
                DrawCycloReducer(inputs)
        except:
            if _ui:
                _ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

class MyCommandExecuteHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()
    def notify(self, args):
        try:
            design = _app.activeProduct
            command = args.firingEvent.sender
            inputs = command.commandInputs
            DrawCycloReducer(inputs)
            saveInputsValues(inputs, USER_CHANGEABLE_ID)

        except:
            if _ui:
                _ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

class MyCommandDestroyHandler(adsk.core.CommandEventHandler):
    def __init__(self):
        super().__init__()
    def notify(self, args):
        try:
            adsk.terminate()
        except:
            if _ui:
                _ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

def run(context):
    global _app, _ui, _unitsMgr
    _ui = None
    try:
        print("addIn start")
        _app = adsk.core.Application.get()
        _ui = _app.userInterface
        _unitsMgr = _app.activeProduct.unitsManager

        cmdDef = _ui.commandDefinitions.itemById(COMMAND_ID)
        if not cmdDef:
            cmdDef = _ui.commandDefinitions.addButtonDefinition(COMMAND_ID, COMMAND_NAME, COMMAND_DESCRIPTION)

        onCommandCreated = MyCommandCreatedHandler()
        cmdDef.commandCreated.add(onCommandCreated)
        handlers.append(onCommandCreated)

        cmdDef.execute()
        adsk.autoTerminate(False)

    except:
        if _ui:
            _ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))