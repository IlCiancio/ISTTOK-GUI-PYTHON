"""
Created on Thu Apr 18 16:09:00 2019

@author: Alessandro Cianciulli, Il_Ciancio
@email: al.cianciulli@studenti.unina.it
"""

import sys
from PyQt5 import QtCore, QtGui, QtWidgets
from GUI_UNINA_5_0 import Ui_MainWindow
from PyQt5.QtCore import Qt

class MainWindow( QtWidgets.QMainWindow ):    
    def __init__(self):
        super(MainWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.Start_Button.clicked.connect(lambda: self.ui.start_on_button_clicked())
        self.ui.Plot_Button.clicked.connect(lambda: self.ui.plot_on_button_clicked())
        self.ui.Plot_Button_Compare.clicked.connect(lambda: self.ui.plot_compare_on_button_clicked())
        self._check_close = True

    def closeEvent(self, event):
        if self._check_close:
            result = QtWidgets.QMessageBox.question(
                self, 'Confirm Close', 'Are you sure you want to close?',
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            if result == QtWidgets.QMessageBox.Yes:
                event.accept()
            elif result == QtWidgets.QMessageBox.No:
                event.ignore()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            self.close()

def main():
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
