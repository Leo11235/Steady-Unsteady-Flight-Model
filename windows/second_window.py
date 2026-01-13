from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QPushButton
from PyQt6.QtCore import Qt

class SecondWindow(QWidget):
    def __init__(self, parent_window):
        super().__init__()
        self.parent_window = parent_window  # Reference to the main window
        self.init_ui()
        
    def init_ui(self):
        # Set window properties
        self.setWindowTitle("Simple App - Window 2")
        self.setGeometry(150, 150, 600, 200)
        
        # Create layout
        layout = QVBoxLayout()
        
        # Create label
        label = QLabel("hi again!")
        label.setStyleSheet("font-size: 24px; padding: 20px;")
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
        # Create back button
        back_button = QPushButton("Go Back")
        back_button.setStyleSheet("padding: 10px; font-size: 16px;")
        back_button.clicked.connect(self.go_back)
        
        # Add widgets to layout
        layout.addWidget(label)
        layout.addWidget(back_button)
        
        # Set layout
        self.setLayout(layout)
        
    def go_back(self):
        """Return to the main window"""
        self.parent_window.show()
        self.close()