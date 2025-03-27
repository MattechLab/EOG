import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QPushButton,
                             QTableWidget, QTableWidgetItem, QComboBox, QLabel)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

class EOG_GUI(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EOG Angle Computation")
        self.setGeometry(100, 100, 800, 600)
        
        # Layout
        self.layout = QVBoxLayout()
        
        # Channel selection dropdowns
        self.channel_x = QComboBox()
        self.channel_y = QComboBox()
        self.channel_x.addItems(["Channel 1", "Channel 2", "Channel 3", "Channel 4"])
        self.channel_y.addItems(["Channel 1", "Channel 2", "Channel 3", "Channel 4"])
        self.layout.addWidget(QLabel("Select Horizontal Channel:"))
        self.layout.addWidget(self.channel_x)
        self.layout.addWidget(QLabel("Select Vertical Channel:"))
        self.layout.addWidget(self.channel_y)
        
        # Compute button
        self.compute_button = QPushButton("Compute Angles")
        self.compute_button.clicked.connect(self.compute_angles)
        self.layout.addWidget(self.compute_button)
        
        # Table to show results
        self.table = QTableWidget(1, 4)
        self.table.setHorizontalHeaderLabels(["Channel X", "Channel Y", "Azimuth Angle", "Elevation Angle"])
        self.layout.addWidget(self.table)
        
        # Matplotlib Canvas
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)
        
        self.setLayout(self.layout)
    
    def compute_angles(self):
        # Placeholder formula for angle computation
        voltage_x = np.random.uniform(-1, 1)
        voltage_y = np.random.uniform(-1, 1)
        azimuth = np.arctan(voltage_x) * (180 / np.pi)
        elevation = np.arctan(voltage_y) * (180 / np.pi)
        
        # Update table
        self.table.setItem(0, 0, QTableWidgetItem(self.channel_x.currentText()))
        self.table.setItem(0, 1, QTableWidgetItem(self.channel_y.currentText()))
        self.table.setItem(0, 2, QTableWidgetItem(f"{azimuth:.2f}"))
        self.table.setItem(0, 3, QTableWidgetItem(f"{elevation:.2f}"))
        
        # Update heatmap
        self.ax.clear()
        heatmap_data = np.random.rand(10, 10)  # Placeholder for actual angle data
        cax = self.ax.imshow(heatmap_data, cmap='jet', interpolation='nearest')
        self.figure.colorbar(cax)
        self.canvas.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = EOG_GUI()
    window.show()
    sys.exit(app.exec())
