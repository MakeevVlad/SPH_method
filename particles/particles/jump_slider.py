from PyQt5 import QtGui
from PyQt5.QtWidgets import QSlider
from PyQt5.QtCore import pyqtSignal


class QJumpSlider(QSlider):
    gaugeClick = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        if hasattr(event, 'type'):
            value = QtGui.QStyle.sliderValueFromPosition(self.minimum(), self.maximum(), event.x(), self.width())
            self.gaugeClick.emit(value)
            self.setValue(value)
