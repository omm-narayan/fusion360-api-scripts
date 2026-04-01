import adsk.core
import adsk.fusion
import traceback

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        design = adsk.fusion.Design.cast(app.activeProduct)

        if not design:
            ui.messageBox('No active Fusion 360 design found.')
            return

        # Warn if design history is on
        if design.designType == adsk.fusion.DesignTypes.ParametricDesignType:
            result = ui.messageBox(
                'Design History is ON.\n\n'
                'Turn it OFF before running for best results:\n'
                'Right-click root component → "Do Not Capture Design History"\n\n'
                'Continue anyway?',
                'Warning',
                adsk.core.MessageBoxButtonTypes.YesNoButtonType
            )
            if result == adsk.core.DialogResults.DialogNo:
                return

        root = design.rootComponent
        moved = 0
        skipped = 0

        # Recursively collect all occurrences
        def get_all_occurrences(occurrences):
            all_occs = []
            for occ in occurrences:
                all_occs.append(occ)
                if occ.childOccurrences.count > 0:
                    all_occs.extend(get_all_occurrences(occ.childOccurrences))
            return all_occs

        all_occs = get_all_occurrences(root.occurrences)

        for i, occ in enumerate(all_occs):
            comp = occ.component

            if comp.bRepBodies.count == 0:
                continue

            for body in comp.bRepBodies:
                try:
                    # copyToComponent places the body in world space automatically
                    # No MoveFeature needed — applying one causes double transform
                    body.copyToComponent(root)
                    moved += 1
                except Exception:
                    skipped += 1

            # Yield control every 5 to keep UI responsive
            if i % 5 == 0:
                adsk.doEvents()

        ui.messageBox(
            f'done!\n\n'
            f'bodies copied to root: {moved}\n'
            f'skipped: {skipped}\n\n'
            f'check positions, then delete the old component folders.'
        )

    except Exception:
        if ui:
            ui.messageBox('Script error:\n' + traceback.format_exc())