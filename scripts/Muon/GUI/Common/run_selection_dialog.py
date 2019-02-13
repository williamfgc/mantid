# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
import PyQt4.QtGui as QtGui

class RunSelectionDialog(QtGui.QDialog):
    def __init__(self, current_runs, parent=None):
        QtGui.QDialog.__init__(self,parent)
        
        layout = QtGui.QVBoxLayout(self)

        current_runs_as_string = [str(run) for run in current_runs]
        self.run_selector_combo = QtGui.QComboBox()
        self.run_selector_combo.addItems(current_runs_as_string)

        layout.addWidget(self.run_selector_combo)

        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtGui.Qt.Horizontal, self)
        
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.setLayout(layout)
    
    def run(self):
        return self.run_selector_combo.getCurrentText()

    @staticmethod
    def get_run(parent = None):
        dialog = RunSelectionDialog(parent)
        result = dialog.exec_()
        run = dialog.run()
        return (run, result == QtGui.QDialog.Accepted)


