import math

import pyqtgraph.opengl as gl
import numpy as np

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from pyqtgraph.opengl import GLViewWidget


class ParticleVisualizer(GLViewWidget):
    posChange = pyqtSignal(int)

    def __init__(self):
        super().__init__()
        self._vis_pos = None
        self._vis_total_frames = None
        self._vis_data = None
        self._vis_scatter = None
        self._vis_n_particles = None
        self._vis_timer = None
        self._vis_is_paused = False
        self.opts['distance'] = 40
        self.setGeometry(0, 110, 1920, 1080)
        self._make_grids()

    def set_data_file(self, path, auto_play=False):
        self._vis_data = np.genfromtxt(path, delimiter=' ')
        self._vis_total_frames = len(self._vis_data)
        if auto_play:
            self.start_playback()

    def start_playback(self):
        if self._vis_is_paused:
            self._vis_timer.start()
            self._vis_is_paused = True
        elif not self._vis_scatter:
            self._init_scatter()
            self._start_animation()

    def pause_playback(self):
        self._vis_timer.stop()
        self._vis_is_paused = True

    def stop_playback(self):
        self.pause_playback()
        self._set_pos(0)

    def randomize_colors(self):
        self._vis_scatter.setData(color=self._get_random_colors())

    def jump_to_pos(self, pos):
        self._set_pos(percent=pos, emit=False)

    def _get_random_colors(self):
        return np.random.rand(self._vis_n_particles, 4) * 2.5

    def _init_scatter(self):
        self._set_pos(frame=0, update=False)
        particles = np.split(self._vis_data[0], len(self._vis_data[0]) / 3)
        self._vis_n_particles = len(particles)

        particles = np.array([tuple(p) for p in particles])

        self._vis_scatter = gl.GLScatterPlotItem(pos=particles, color=self._get_random_colors(), size=10)
        self.addItem(self._vis_scatter)

    def _make_grids(self):
        gx = gl.GLGridItem()
        gx.rotate(90, 0, 1, 0)
        gx.translate(-10, 0, 0)
        self.addItem(gx)

        gy = gl.GLGridItem()
        gy.rotate(90, 1, 0, 0)
        gy.translate(0, -10, 0)
        self.addItem(gy)

        gz = gl.GLGridItem()
        gz.translate(0, 0, -10)
        self.addItem(gz)

    def _update_plots(self):
        frame_data = self._vis_data[self._vis_pos]
        particles = np.split(frame_data, len(frame_data) / 3)
        particles = np.array([tuple(p) for p in particles])
        self._vis_scatter.setData(pos=particles)

    def _set_pos(self, frame=None, percent=None, emit=True, update=True):
        if percent is not None:
            frame = percent * int(self._vis_total_frames / 100)
        self._vis_pos = frame
        if update:
            self._update_plots()
        if emit:
            self.posChange.emit(int(100 * self._vis_pos / self._vis_total_frames))  # percent

    def _increment_pos(self):
        self._set_pos(frame=self._vis_pos + 1 if self._vis_pos < len(self._vis_data) - 1 else 0)

    def _start_animation(self):
        self._vis_timer = QtCore.QTimer()
        self._vis_timer.timeout.connect(self._increment_pos)
        self._vis_timer.start(0.05)
