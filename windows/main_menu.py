from PyQt6.QtWidgets import QMainWindow, QMessageBox
from PyQt6 import uic
import os

class MainMenu(QMainWindow):
    def __init__(self):
        super().__init__()
        
        # Load the .ui file
        current_dir = os.path.dirname(__file__)
        ui_file = os.path.join(current_dir, "..", "ui", "main_menu.ui")
        uic.loadUi(ui_file, self)
        
        # Connect buttons to functions
        self.connect_buttons()
        
    def connect_buttons(self):
        """Connect all button signals to their slots"""
        self.pushButton_1.clicked.connect(self.on_button1_clicked)
        self.pushButton_2.clicked.connect(self.on_button2_clicked)
        self.pushButton_3.clicked.connect(self.on_button3_clicked)
        self.pushButton_4.clicked.connect(self.on_button4_clicked)
        
    def on_button1_clicked(self):
        QMessageBox.information(self, "Button 1", "Button 1 was clicked!")
        # Add your logic here
        
    def on_button2_clicked(self):
        QMessageBox.information(self, "Button 2", "Button 2 was clicked!")
        # Add your logic here
        
    def on_button3_clicked(self):
        QMessageBox.information(self, "Button 3", "Button 3 was clicked!")
        # Add your logic here
        
    def on_button4_clicked(self):
        QMessageBox.information(self, "Button 4", "Button 4 was clicked!")
        # Add your logic here