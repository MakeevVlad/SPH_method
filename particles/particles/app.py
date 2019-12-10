import os
import sys
from pathlib import Path

from PyQt5 import Qt, QtCore
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog, QAction, QSlider, QVBoxLayout, QPushButton, \
    QGridLayout, QWidget, QHBoxLayout, QToolButton, QStatusBar

from particles.jump_slider import QJumpSlider
from particles.vis import ParticleVisualizer
from settings import _resource_path, SAMPLE_DATA_FILE, BASE_DIR, MAX_RECENTS, RECENTS_PATH, MAX_RECENT_LENGTH


class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)
        self._set_title()
        self.resize(1024, 768)
        self._toolbar = None
        self._recents_menu = None

        self._init_menu()

        self._init_toolbar()
        self._init_widgets()

    def _init_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu('Файл')
        file_menu.aboutToShow.connect(self._check_recents_enabled)

        open_action = QAction('Открыть...', self)
        open_action.triggered.connect(self._open_file_dialog)
        file_menu.addAction(open_action)

        self._recents_menu = file_menu.addMenu('Открыть недавние')
        self._recents_menu.aboutToShow.connect(self._set_recents_menu)

        load_sample = QAction('Загрузить демо', self)
        load_sample.triggered.connect(self._load_sample)
        file_menu.addAction(load_sample)

        close_action = QAction('Выход', self)
        close_action.triggered.connect(self.close)
        file_menu.addAction(close_action)

    def _init_toolbar(self):
        toolbar = self.addToolBar('Управление')
        toolbar.setEnabled(False)

        play_action = QAction(QIcon(_resource_path('icon-play.svg')), 'Воспроизведение', self)
        play_action.triggered.connect(self._start_playback)
        toolbar.addAction(play_action)

        pause_action = QAction(QIcon(_resource_path('icon-pause.svg')), 'Пауза', self)
        pause_action.triggered.connect(self._pause_playback)
        toolbar.addAction(pause_action)

        stop_action = QAction(QIcon(_resource_path('icon-stop.svg')), 'Стоп', self)
        stop_action.triggered.connect(self._stop_playback)
        toolbar.addAction(stop_action)

        recolor_action = QAction(QIcon(_resource_path('icon-recolor.svg')), 'recolor', self)
        recolor_action.triggered.connect(self._recolor)
        toolbar.addAction(recolor_action)

        self._toolbar = toolbar

    def _init_widgets(self):
        widget = QWidget(self)
        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(0, 0, 0, 0)

        self._vis = ParticleVisualizer()
        self._vis.posChange.connect(self._on_pos_change)
        vbox.addWidget(self._vis, stretch=1)

        hbox = QHBoxLayout(self)
        slider = self._slider = QJumpSlider(QtCore.Qt.Horizontal)
        slider.setTickPosition(QJumpSlider.TicksBothSides)
        slider.setMaximum(0)
        slider.setMaximum(99)
        slider.setTickInterval(1)
        slider.setEnabled(False)
        slider.sliderMoved.connect(self._jump_to_pos)
        slider.gaugeClick.connect(self._jump_to_pos)
        hbox.addWidget(self._slider)
        hbox.setContentsMargins(10, 5, 10, 10)
        vbox.addLayout(hbox)

        widget.setLayout(vbox)
        self.setCentralWidget(widget)

    def _get_recents(self):
        if not os.path.isfile(RECENTS_PATH):
            Path(RECENTS_PATH).touch()
        recents = []
        with open(RECENTS_PATH, 'r', encoding='utf-8') as file:
            for recent_path in file.readlines():
                recent_path = recent_path.strip()
                recents.append((recent_path, os.path.isfile(recent_path)))
        return recents

    def _add_path_to_recents(self, path):
        recents = [rec[0] for rec in self._get_recents() if rec[0] != path]
        recents = [path, *recents[:MAX_RECENTS - 1]]
        with open(RECENTS_PATH, 'w', encoding='utf-8') as file:
            for recent in recents:
                file.write(f'{recent}\n')

    def _check_recents_enabled(self):
        self._recents_menu.setEnabled(len(self._get_recents()))

    def _set_recents_menu(self):
        self._recents_menu.clear()
        recents = self._get_recents()
        for recent_path, is_valid in recents:
            name = recent_path if len(recent_path) <= MAX_RECENT_LENGTH else f'...{recent_path[-MAX_RECENT_LENGTH:]}'
            action = QAction(name, self)
            action.setEnabled(is_valid)
            action.triggered.connect(lambda _, path=recent_path: self._load_data_file(path))
            self._recents_menu.addAction(action)

    def _open_file_dialog(self):
        filename = QFileDialog.getOpenFileName()
        path = filename[0]
        if path:
            self._load_data_file(path)
            self._add_path_to_recents(path)

    def _load_sample(self):
        self._load_data_file(SAMPLE_DATA_FILE)

    def _recolor(self):
        self._vis.randomize_colors()

    def _start_playback(self):
        self._vis.start_playback()

    def _pause_playback(self):
        self._vis.pause_playback()

    def _stop_playback(self):
        self._vis.stop_playback()

    def _load_data_file(self, path):
        self._set_title(os.path.basename(path))
        self._vis.set_data_file(path, auto_play=True)
        self._toolbar.setEnabled(True)
        self._slider.setEnabled(True)

    def _set_title(self, filename=None):
        self.setWindowTitle('Particles' if not filename else f'Particles [{filename}]')

    def _jump_to_pos(self, pos):
        self._vis.jump_to_pos(pos)

    def _on_pos_change(self, pos):
        self._slider.setValue(pos)


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    app.exec_()
