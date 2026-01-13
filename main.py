import os
import sys
from PyQt6.QtWidgets import QApplication
from windows.main_menu import MainMenu

def main():
    app = QApplication(sys.argv)
    
    # Create and show the main menu
    window = MainMenu()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()