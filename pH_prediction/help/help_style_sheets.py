# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:09:36 2019

@author: jeromlu2
"""

from PyQt5.QtWidgets import QDialog, QApplication,QVBoxLayout, QPushButton

class Test(QDialog):
    """A sample widget to show dyanmic properties with stylesheets"""
    def __init__(self):
        # Do the usual
        super(Test, self).__init__()
        self.setWindowTitle('Style Sheet Test')
        layout = QVBoxLayout(self)

        # Set the stylesheet
        self.setStyleSheet("""
            QPushButton[Test=true] {
                border: 2px solid #8f8f91;
                border-radius: 6px;
                background-color: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                                      stop: 0 #f6f7fa, stop: 1 #dadbde);
                min-width: 80px;
            }

            QPushButton#StyledButton[Test=true] {
                color: #F00;
                background-color: #000;
            }
                           """)

        # Create the button
        btn = QPushButton('Click Me')
        btn.setProperty('Test', True)
        btn.setObjectName('StyledButton')
        btn.clicked.connect(lambda: self.toggle(btn))
        layout.addWidget(btn)

        btn2 = QPushButton('Click Me')
        btn2.setProperty('Test', True)
        btn2.clicked.connect(lambda: self.toggle(btn2))
        layout.addWidget(btn2)

    def toggle(self, widget):
        # Query the attribute
        isTest = widget.property('Test')
        widget.setProperty('Test', not isTest)

        # Update the style
        widget.setStyle(widget.style())

if __name__ == '__main__':
    app = QApplication([])
    dlg = Test()
    dlg.show()
    app.exec_()