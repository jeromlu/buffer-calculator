# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:29:30 2019

@author: jeromlu2
"""

import sys

from PyQt5.QtWidgets import QApplication, QDialog

from PyQt5.QtWidgets import QLabel, QLineEdit, QBoxLayout
from PyQt5.QtWidgets import QWidget, QComboBox

from PyQt5.QtGui import QDoubleValidator, QValidator


#Global constants
LEFT, ABOVE = range(2)

class LabelledLineEdit(QWidget):
    
    def __init__(self, labelText = '', position=LEFT,
        parent=None):
        super(LabelledLineEdit, self).__init__(parent)
        self.label = QLabel(labelText)
        self.line_edit = QLineEdit()
        self.label.setBuddy(self.line_edit)
        layout = QBoxLayout(QBoxLayout.LeftToRight \
        if position == LEFT else QBoxLayout.TopToBottom)
        layout.addWidget(self.label)
        layout.addWidget(self.line_edit)
        self.setLayout(layout)
        
        
class LabelledComboBox(QWidget):
    
    def __init__(self, labelText = '', position=LEFT,
        parent=None):
        super(LabelledComboBox, self).__init__(parent)
        self.label = QLabel(labelText)
        self.combo_box = QComboBox()
        self.label.setBuddy(self.combo_box)
        layout = QBoxLayout(QBoxLayout.LeftToRight \
        if position == LEFT else QBoxLayout.TopToBottom)
        layout.addWidget(self.label)
        layout.addWidget(self.combo_box)
        self.setLayout(layout)
        
class DoubleLLE(LabelledLineEdit):
    
    def __init__(self, labelText = '', position=LEFT, parent=None):
        super(DoubleLLE, self).__init__(labelText, position, parent)    
    
        validator = QDoubleValidator()
        self.line_edit.setValidator(validator)
        self.line_edit.textChanged.connect(self.check_state)
        self.line_edit.textChanged.emit(self.line_edit.text())
        
    def check_state(self, *args, **kwargs):
        sender = self.sender()
        validator = sender.validator()
        state = validator.validate(sender.text(), 0)[0]
        if state == QValidator.Acceptable:
            color = '#c4df9b' # green
        elif state == QValidator.Intermediate:
            color = '#fff79a' # yellow
        else:
            color = '#f6989d' # red
        sender.setStyleSheet('QLineEdit { background-color: %s }' % color)
        
class TestDialog(QDialog):
    
    def __init__(self, parent = None):
        super(TestDialog, self).__init__(parent)
        
        dll = DoubleLLE('test', parent = self)
        
    def closeEvent(self, evnt):
        QApplication.quit()
        
        
if __name__ == "__main__":
    
    app = QApplication(sys.argv)
    
    test_form = TestDialog()
    test_form.show()
    app.exec_()
    