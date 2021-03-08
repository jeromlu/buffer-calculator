# -*- coding: utf-8 -*-
#
# Created on Tue Mar 05 2019
#
# Copyright (c) 2021 Your Company
# Name: Luka Jeromel
#
# ******************************Python imports***********************************
import sys
import time
import io
import traceback
import os

# ******************************PyQt5 imports*************************************
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QAction
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem, QVBoxLayout
from PyQt5.QtWidgets import QHeaderView, QPushButton, QLabel
from PyQt5.QtWidgets import QHBoxLayout, QTabWidget, QFrame
from PyQt5.QtWidgets import QAbstractScrollArea, QDialog
from PyQt5.QtWidgets import QMessageBox, QSizePolicy, QCheckBox

from PyQt5.QtGui import QIcon, QTextDocument

from PyQt5.QtCore import Qt

from PyQt5.QtPrintSupport import QPrintDialog, QPrinter, QPrintPreviewDialog

# ******************************Other third party imports************************
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

# ******************************My modules***************************************
import custom_widgets as cw
import pH_calculation_Rushd
import output_string
import qrc_resources  # noqa F401

# exe location
# determine if application is a script file or frozen exe
if getattr(sys, "frozen", False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)


# Global constants
LEFT, ABOVE = range(2)
ACID_NAMES = [
    "Acetate",
    "Phosphate",
    "Succinate",
    "Tris",
    "Citrate",
    "Malate",
    "Arginine",
    "Histidine",
    "Cysteine",
    "EDTA",
    "NaCl",
]

CHEMICALS = [
    "H2O",
    "Acetic acid",
    "Sodium acetate, anhydrous",
    "o-Phosphoric acid >85%",
    "Sodium di-hydrogen phosphate monohydrate (NaH2PO4)",
    "di-Sodium hydrogen phosphate (Na2HPO4), anhydrous",
    "tri-Sodium phosphate (Na3HPO4), anhydrous",
    "succinic acid",
    "Tris(hydroxymethyl)aminomethane",
    "Tris(hydroxymethyl)aminomethan hydrochloride",
    "Citric acid monohydrate",
    "Malic acid, anhydrous",
    "L-Arginine monohydrochloride",
    "L-Histidine",
    "L-Histidine monohydrochloride monohydrate",
    "L-cysteine hydrochloride monohydrate",
    "EDTA-disodium salt (Titriplex III) di hydrate",
    "Sodium chloride",
    "sodium hydroxide solution 10.0 M",
    "Hydrochloric acid 25%",
]

PHOSPHATE = {"NaH2PO4": 1, "Na2HPO4": 2, "Na3HPO4": 3}

COLOR_PH = "tab:blue"
COLOR_COND = "tab:green"


class BufferCalculatorUI(QMainWindow):
    def __init__(self, parent=None):
        super(BufferCalculatorUI, self).__init__(parent)

        # list containing the mass of each component
        self.composition = []
        # dictionary of lists [temperature, ionic strength, pH, conductivity]
        self.temp_dep = {"Temp": [], "IS": [], "pH": [], "Cond": []}
        self.plot_lns = [None, None]
        self.print_document = None
        self.printer = None

        # Create the UI
        self.create_main_window()
        self.create_menu_bar()

        # connections
        self.pb_calc_composition.clicked.connect(self.calculate_buff_composition)
        self.lle_buffer_name.line_edit.textChanged.connect(self.update_buffer_title)
        # self.pb_export.clicked.connect(self.mpl_canvas.draw)
        self.pb_preview.clicked.connect(self.preview_recepie)
        self.pb_print.clicked.connect(self.print_recepie)
        self.table_molarities.cellChanged.connect(self.create_name)

    def create_main_window(self):

        self.main_frame = QWidget(self)

        # Table widget - could create separate class...
        stylesheet = "alternate-background-color: rgb(212, 247, 226)"
        self.table_molarities = QTableWidget()
        self.table_molarities.setStyleSheet(stylesheet)
        self.table_molarities.setAlternatingRowColors(True)
        # set row count
        self.table_molarities.setRowCount(len(ACID_NAMES))
        # set column count
        self.table_molarities.setColumnCount(2)
        self.table_molarities.setHorizontalHeaderLabels(["Acid", "Amount [mM]"])
        self.table_molarities.horizontalHeaderItem(1).setTextAlignment(Qt.AlignHCenter)
        self.table_molarities.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
        self.table_molarities.setSortingEnabled(False)
        self.table_molarities.verticalHeader().setVisible(False)
        self.table_molarities.horizontalHeader().setSectionsMovable(False)
        self.table_molarities.horizontalHeader().setSectionResizeMode(QHeaderView.Fixed)
        self.table_molarities.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        for i, name in enumerate(ACID_NAMES):
            item = QTableWidgetItem(name)
            item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            self.table_molarities.setItem(i, 0, item)
            item = QTableWidgetItem()
            if i == 1:
                item.setData(Qt.EditRole, 50.0)
            else:
                item.setData(Qt.EditRole, 0.0)
            self.table_molarities.setItem(i, 1, item)

        lle_width = 100
        self.lle_volume = cw.DoubleLLE("Volume [L]", ABOVE)
        self.lle_volume.setFixedWidth(lle_width)
        self.lle_volume.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.lle_volume.line_edit.insert("1")
        self.lle_pH = cw.DoubleLLE("Buffer pH", ABOVE)
        self.lle_pH.setFixedWidth(lle_width)
        self.lle_pH.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.lle_pH.line_edit.insert("7")
        self.lcb_phosphate = cw.LabelledComboBox("Phosphate chemical", ABOVE)
        self.lcb_phosphate.combo_box.insertItems(0, PHOSPHATE.keys())
        self.lcb_phosphate.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Expanding)
        self.lle_buffer_name = cw.LabelledLineEdit("Name of a buffer", ABOVE)
        self.lle_buffer_name.setFixedWidth(300)

        self.cb_auto_name = QCheckBox("automatic naming")

        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tab_composition = QWidget()
        self.tab_temperature = QWidget()
        self.tabs.resize(200, 500)

        # Add tabs
        self.tabs.addTab(self.tab_composition, "Buffer composition")
        self.tabs.addTab(self.tab_temperature, "Temperature dependence")

        # Create compostion tab
        self.pb_export = QPushButton("E&xport composition")
        self.pb_preview = QPushButton("Print previe&w")
        self.pb_print = QPushButton("&Print")
        hbox_export = QHBoxLayout()
        hbox_export.addWidget(self.pb_export)
        hbox_export.addWidget(self.pb_preview)
        hbox_export.addWidget(self.pb_print)
        hbox_export.addStretch(1)

        # Table widget - could create separate class...
        stylesheet = "alternate-background-color: rgb(210, 255, 255)"
        self.table_chemicals = QTableWidget()
        self.table_chemicals.setStyleSheet(stylesheet)
        self.table_chemicals.setAlternatingRowColors(True)
        # set row count
        self.table_chemicals.setRowCount(len(CHEMICALS))
        # set column count
        self.table_chemicals.setColumnCount(2)
        self.table_chemicals.setHorizontalHeaderLabels(["Amount [g]", "Chemical"])
        self.table_chemicals.horizontalHeaderItem(1).setTextAlignment(Qt.AlignHCenter)
        self.table_chemicals.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
        self.table_chemicals.setSortingEnabled(False)
        self.table_chemicals.verticalHeader().setVisible(False)
        self.table_chemicals.horizontalHeader().setSectionsMovable(False)
        self.table_chemicals.horizontalHeader().setSectionResizeMode(QHeaderView.Fixed)
        for i, name in enumerate(CHEMICALS):
            item = QTableWidgetItem(name)
            item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            self.table_chemicals.setItem(i, 1, item)
            item = QTableWidgetItem()
            item.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            self.table_chemicals.setItem(i, 0, item)
        self.table_chemicals.resizeColumnsToContents()

        # additional info about the
        style_sheet = "QLabel { background-color : lightgrey; color : black; }"
        cond_label = QLabel("Conductivity")
        self.cond_label = QLabel("Conductivity")
        self.cond_label.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        self.cond_label.setMinimumHeight(30)
        self.cond_label.setStyleSheet(style_sheet)
        IS_label = QLabel("Ionic strength")
        self.IS_label = QLabel("Ionic strength")
        self.IS_label.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        self.IS_label.setMinimumHeight(30)
        self.IS_label.setStyleSheet(style_sheet)

        vbox_add_info = QVBoxLayout()
        vbox_add_info.addWidget(cond_label)
        vbox_add_info.addWidget(self.cond_label)
        vbox_add_info.addWidget(IS_label)
        vbox_add_info.addWidget(self.IS_label)
        vbox_add_info.addStretch(1)

        hbox_chemicals = QHBoxLayout()
        hbox_chemicals.addWidget(self.table_chemicals)
        hbox_chemicals.addLayout(vbox_add_info)
        hbox_chemicals.addStretch(1)

        self.tab_composition.layout = QVBoxLayout()
        self.tab_composition.layout.addLayout(hbox_chemicals)
        self.tab_composition.layout.addLayout(hbox_export)
        self.tab_composition.setLayout(self.tab_composition.layout)

        # create figure tab (temperature dependence)
        self.figure = Figure(figsize=(400, 300), dpi=100)

        # create figure widget with navigation toolbar
        self.mpl_canvas = FigureCanvas(self.figure)
        self.mpl_canvas.setParent(self.tab_temperature)
        self.navigation_bar = NavigationToolbar(
            self.mpl_canvas, self.tab_temperature, coordinates=True
        )

        # create axes to plot on, and initialize some parameters
        # pH axis
        self.ax_pH = self.mpl_canvas.figure.add_subplot(111)
        self.ax_pH.set_xlim(14, 36)
        self.ax_pH.set_xlabel("Temperature [C]")
        self.ax_pH.set_ylabel("pH", color=COLOR_PH)
        self.ax_pH.tick_params(axis="y", labelcolor=COLOR_PH)

        # cond axis
        self.ax_cond = self.ax_pH.twinx()
        self.ax_cond.set_ylabel("Conductivity [mS/cm]", color=COLOR_COND)
        self.ax_cond.tick_params(axis="y", labelcolor=COLOR_COND)
        # vbox_temp = QVBoxLayout()
        # vbox_temp.addWidget(self.mpl_canvas)
        # vbox_temp.addWidget(self.navigation_bar)

        self.tab_temperature.layout = QVBoxLayout()
        self.tab_temperature.layout.addWidget(self.mpl_canvas)
        self.tab_temperature.layout.addWidget(self.navigation_bar)
        self.tab_temperature.setLayout(self.tab_temperature.layout)

        self.tabs.setCurrentIndex(1)

        self.pb_calc_composition = QPushButton("Ca&lculate\ncomposition")

        vbox_lle = QVBoxLayout()
        vbox_lle.setAlignment(Qt.AlignLeft)
        vbox_lle.addWidget(self.lle_pH)
        vbox_lle.addWidget(self.lle_volume)
        vbox_lle.addWidget(self.lcb_phosphate)
        vbox_lle.addWidget(self.pb_calc_composition)
        vbox_lle.addStretch(1)

        vbox_table = QVBoxLayout()
        vbox_table.setAlignment(Qt.AlignCenter)
        vbox_table.addWidget(self.table_molarities)
        vbox_table.addStretch(1)

        hbox2 = QHBoxLayout()
        hbox2.addLayout(vbox_table)
        hbox2.addLayout(vbox_lle)

        vbox_left = QVBoxLayout()
        vbox_left.addWidget(self.lle_buffer_name)
        vbox_left.addLayout(hbox2)
        vbox_left.addStretch(1)
        vbox_left.addWidget(self.cb_auto_name)

        hbox_both = QHBoxLayout()
        hbox_both.addLayout(vbox_left)
        hbox_both.addWidget(self.tabs)

        # vbox_all = QVBoxLayout()
        # vbox_all.addLayout(hbox_table)
        # vbox_all.addLayout(vbox_left)

        self.setCentralWidget(self.main_frame)
        self.main_frame.setLayout(hbox_both)

    def create_menu_bar(self):

        # file menu
        self.file_menu = self.menuBar().addMenu("&File")

        file_quit_action = self.create_action(
            "&Close app",
            shortcut="Ctrl+Q",
            slot=QApplication.quit,
            tip="Close the application",
            icon=None,
        )

        self.add_actions(self.file_menu, [file_quit_action])

    def calculate_buff_composition(self):
        out_array = []
        for i in range(len(ACID_NAMES)):
            item = self.table_molarities.item(i, 1)
            out_array.append(item.text())
        out_array.append(self.lle_pH.line_edit.text())
        key = self.lcb_phosphate.combo_box.currentText()
        out_array.append(str(PHOSPHATE[key]))
        out_array.append(self.lle_volume.line_edit.text())
        if out_array[-1] == "":
            QMessageBox.warning(
                self,
                "Volume error",
                "If you are preparing buffer of no volume,"
                "then you are already finished.\n\nCongratiolations!! :)",
            )
            return
        communicator = pH_calculation_Rushd.FortranCommunication(
            folder=application_path + "/testing/"
        )
        communicator.write_parameters(out_array)
        communicator.run_subprocess()
        if not communicator.read_data(self):
            QMessageBox.critical(
                self,
                "Error",
                "Check your values,\n"
                "something was wrong with calculation,\n"
                "it could be too high molarity.",
            )
            self.clear_data()
            return
        self.update_composition()
        self.update_temp_dep()

    def create_name(self):

        if not self.cb_auto_name.isChecked():
            return

        # maximal number of characters
        max_char = 250
        buff_name = ""
        # iterate thorugh the all possible components
        for i, a_name in enumerate(ACID_NAMES):
            item = self.table_molarities.item(i, 1)
            if float(item.text()) > 0:
                buff_name += item.text()
                buff_name = buff_name + " mM "
                buff_name += a_name
                buff_name += " "
        buff_name = buff_name.rstrip()
        self.lle_buffer_name.line_edit.setText(buff_name[:max_char])

    def update_composition(self):
        for i, value in enumerate(self.composition):
            item = QTableWidgetItem(value)
            item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.table_chemicals.setItem(i, 0, item)

    def update_temp_dep(self):

        try:
            x = self.temp_dep["Temp"]
            y_pH = self.temp_dep["pH"]
            y_cond = self.temp_dep["Cond"]
            self.cond_label.setText("{0:.2f} mS/cm".format(y_cond[6]))
            self.IS_label.setText("{0:.2f} g/L".format(self.temp_dep["IS"][6]))

            target_pH = float(self.lle_pH.line_edit.text())

            self.ax_pH.set_title(self.lle_buffer_name.line_edit.text())

            if self.plot_lns[0]:
                self.plot_lns[0].set_data(x, y_pH)
                self.plot_lns[1].set_data(x, y_cond)
                self.ax_pH.set_ylim(target_pH - 2, target_pH + 2)
                self.ax_cond.set_ylim(min(y_cond), max(y_cond))

            else:
                self.ax_pH.set_ylim(target_pH - 2, target_pH + 2)
                self.ax_cond.set_ylim(min(y_cond), max(y_cond))
                (self.plot_lns[0],) = self.ax_pH.plot(
                    x, y_pH, label="Buffer pH", color=COLOR_PH, linewidth=2
                )

                (self.plot_lns[1],) = self.ax_cond.plot(
                    x,
                    y_cond,
                    label="Buffer conductivity",
                    color=COLOR_COND,
                    linewidth=2,
                )

            self.mpl_canvas.draw()

            # self.mpl_canvas.flush_events()

        except Exception as e:
            print(e)

    def update_buffer_title(self):
        self.ax_pH.set_title(self.lle_buffer_name.line_edit.text())
        self.mpl_canvas.draw()

    def clear_data(self):

        # list containong the mass of each component
        self.composition = []
        # dictionary of lists [temperature, ionic strength, pH, conductivity]
        self.temp_dep = {"Temp": [], "IS": [], "pH": [], "Cond": []}

    def create_print_document(self):
        name = self.lle_buffer_name.line_edit.text()
        if name == "":
            name = "Unnamed buffer"
        pH = self.lle_pH.line_edit.text()
        cond = self.cond_label.text()
        text = output_string.create_string(name, CHEMICALS, self.composition, pH, cond)
        self.print_document = QTextDocument()
        self.print_document.setHtml(text)

    def print_recepie(self):

        self.create_print_document()

        if self.print_document is None:
            return
        if self.printer is None:
            self.printer = QPrinter(QPrinter.HighResolution)
            self.printer.setPageSize(QPrinter.Letter)
        dialog = QPrintDialog()
        if dialog.exec_() == QDialog.Accepted:
            self.print_document.print_(dialog.printer())

    def preview_recepie(self):
        self.create_print_document()
        dialog = QPrintPreviewDialog()
        dialog.paintRequested.connect(self.print_document.print_)
        dialog.exec_()

    def closeEvent(self, evnt):
        QApplication.quit()

    # ********************************************* HELPER FUNCTIONS

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(
        self,
        text,
        slot=None,
        shortcut=None,
        icon=None,
        tip=None,
        checkable=False,
        signal="triggered()",
    ):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon("./icons/{}.png".format(icon)))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action

    def print_err(self):
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        err_msg = "{0}:\n{1}\nError occurred in file: {2}".format(
            exc_type, exc_obj, fname
        )
        # QMessageBox.critical(self, 'Error - see below', err_msg)


def excepthook(excType, excValue, tracebackobj):
    """Global function to catch unhandled exceptions.

    @param excType exception type
    @param excValue exception value
    @param tracebackobj traceback object
    """
    separator = "-" * 80
    logFile = "simple.log"
    notice = (
        """An unhandled exception occurred. Please report the problem\n"""
        """using the error reporting dialog or via email to <{0}>.\n"""
        """A log has been written to {1}.\n\nError information:\n""".format(
            "luka.jeromel@novartis.com", ""
        )
    )
    versionInfo = "0.0.1"
    timeString = time.strftime("%Y-%m-%d, %H:%M:%S")

    tbinfofile = io.StringIO()
    traceback.print_tb(tracebackobj, None, tbinfofile)
    tbinfofile.seek(0)
    tbinfo = tbinfofile.read()
    errmsg = "%s: \n%s" % (str(excType), str(excValue))
    sections = [separator, timeString, separator, errmsg, separator, tbinfo]
    msg = "\n".join(sections)
    try:
        f = open(logFile, "w")
        f.write(msg)
        f.write(versionInfo)
        f.close()
    except IOError:
        pass
    errorbox = QMessageBox()
    errorbox.setText(str(notice) + str(msg) + str(versionInfo))
    errorbox.exec_()


def main():
    if not QApplication.instance():
        app = QApplication(sys.argv)
    else:
        app = QApplication.instance()
    mainWin = BufferCalculatorUI()
    app.setWindowIcon(QIcon(":/main_window_icon.png"))
    mainWin.setWindowTitle("Buffer maker")
    mainWin.resize(1000, 600)
    mainWin.show()
    app.exec_()


if __name__ == "__main__":

    sys.excepthook = excepthook
    main()
