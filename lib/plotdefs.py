"""
Handy matplotlib.pyplot settings.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pp
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.colors as cm
import matplotlib.font_manager as fm
from cycler import cycler
import copy

color_cycles = {
    "matlab": [
        "#0072bd",
        "#d95319",
        "#edb120",
        "#7e2f8e",
        "#77ac30",
        "#4dbeee",
        "#a2142f",
    ],
    "python": [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ],
    "best": [
        "#1f77b4",  # blue
        "#d95319",  # auburn
        "#edb120",  # canary
        "#7e2f8e",  # purple
        "#46add9",  # cyan
        "#ff7f0e",  # tangerine
        "#3d786e",  # dark seafoam
        "#505050",  # gray
        "#a2142f",  # burgundy
        "#bf7878",  # dark rose
    ]
}

rcdefs = {
    "axes.grid"             : True,
    "axes.grid.which"       : "both",
    "axes.linewidth"        : 0.8,
    "axes.prop_cycle"       : cycler(color=color_cycles["best"]),
    "axes.titlesize"        : "medium",
    "errorbar.capsize"      : 1.0,
    "figure.dpi"            : 500.0,
    "figure.figsize"        : [3.375, 2.100],
    "font.size"             : 10,
    "grid.color"            : "#d8d8d8",
    "grid.linewidth"        : 0.5,
    "image.composite_image" : False,
    "image.cmap"            : "bone",
    "legend.borderaxespad"  : 0.25,
    "legend.borderpad"      : 0.3,
    "legend.fancybox"       : False,
    "legend.fontsize"       : 9, # "small",
    "legend.framealpha"     : 0.8,
    "legend.handlelength"   : 1.2,
    "legend.handletextpad"  : 0.4,
    "legend.labelspacing"   : 0.25,
    "lines.linewidth"       : 1.0,
    "lines.markeredgewidth" : 1.0,
    "lines.markersize"      : 2.5,
    "savefig.bbox"          : "tight",
    "savefig.pad_inches"    : 0.01,
    "text.latex.preamble"   : r"\usepackage{physics}\usepackage{siunitx}\usepackage{amsmath}",
    "xtick.direction"       : "in",
    "xtick.major.size"      : 2.0,
    "xtick.minor.size"      : 1.5,
    "ytick.direction"       : "in",
    "ytick.major.size"      : 2.0,
    "ytick.minor.size"      : 1.5,
}
for key in rcdefs:
    pp.rcParams[key] = rcdefs[key]

hot_cold_colors = [
    (0.000, "#101010"),
    (0.100, "#3f119d"),
    (0.350, "#3967d0"),
    (0.500, "#f0f0f0"),
    (0.625, "#f1b931"),
    (1.000, "#dd0000"),
]
hot_cold = cm.LinearSegmentedColormap.from_list("hot-cold", hot_cold_colors)

plasma_colors = [
    (0.000, "#000000"),
    (0.450, "#3b4568"),
    (0.600, "#586186"),
    (0.700, "#939cc4"),
    (1.000, "#ffffff"),
]
plasma = cm.LinearSegmentedColormap.from_list("plasma", plasma_colors)

cyborg_colors = [
    (0.000, "#101010"),
    (0.100, "#3967d0"),
    (1.000, "#dd0000"),
]
cyborg = cm.LinearSegmentedColormap.from_list("cyborg", cyborg_colors)

vibrant_colors = [
    (0.000, "#101010"),
    (0.050, "#012d5e"),
    (0.125, "#0039a7"),
    (0.250, "#1647cf"),
    (0.375, "#6646ff"),
    (0.500, "#bc27ff"),
    (0.600, "#dc47af"),
    (0.800, "#f57548"),
    (0.900, "#f19e00"),
    (0.950, "#fbb800"),
    (1.000, "#fec800"),
]
vibrant = cm.LinearSegmentedColormap.from_list("vibrant", vibrant_colors)

artsy_colors = [
    (0.000, "#1f0109"),
    (0.034, "#1f0110"),
    (0.069, "#230211"),
    (0.103, "#250816"),
    (0.138, "#270b1b"),
    (0.172, "#250f1d"),
    (0.207, "#251521"),
    (0.241, "#251a25"),
    (0.276, "#2c1b28"),
    (0.310, "#271d2b"),
    (0.345, "#24202d"),
    (0.379, "#232632"),
    (0.414, "#212d32"),
    (0.448, "#1e343c"),
    (0.483, "#173e44"),
    (0.517, "#17464a"),
    (0.552, "#104a49"),
    (0.586, "#0e5553"),
    (0.621, "#00635f"),
    (0.655, "#007065"),
    (0.690, "#007a6d"),
    (0.724, "#0e8476"),
    (0.759, "#1c8c7d"),
    (0.793, "#219581"),
    (0.828, "#2f9f8a"),
    (0.862, "#49a890"),
    (0.897, "#60b89d"),
    (0.931, "#7ec8a9"),
    (0.966, "#9ad6b4"),
    (1.000, "#bce6bf"),
]
artsy = cm.LinearSegmentedColormap.from_list("artsy", artsy_colors)

pix_colors = [
    (0.000, "#0d2b45"),
    (0.143, "#16334d"),
    (0.286, "#544e68"),
    (0.429, "#8d697a"),
    (0.571, "#d08159"),
    (0.714, "#ffaa5e"),
    (0.857, "#ffd4a3"),
    (1.000, "#ffecd6"),
]
pix = cm.LinearSegmentedColormap.from_list("pix", pix_colors)

def figure3D(*fig_args, **fig_kwargs):
    fig = pp.figure(*fig_args, **fig_kwargs)
    ax = p3.Axes3D(fig)
    return fig, ax

def set_font(path, name):
    fe = fm.FontEntry(fname=path, name=name)
    fm.fontManager.ttflist.insert(0, fe)
    pp.rcParams["font.family"] = fe.name

def use_tex(u=True):
    pp.rcParams["text.usetex"] = u
    if u:
        pp.rcParams["font.serif"] = ["Computer Modern Roman"]
        pp.rcParams["font.sans-serif"] = ["Computer Modern Sans Serif"]
        pp.rcParams["font.family"] = ["serif"]
    else:
        pp.rcParams["font.serif"] = ["DejaVu Serif"]
        pp.rcParams["font.sans-serif"] = ["DejaVu Sans"]
        pp.rcParams["font.family"] = ["sans-serif"]

def grid(onoff=True, axes=None):
    if axes:
        axes.minorticks_on()
        if onoff:
            axes.grid(onoff, "major", color="#d8d8d8")
            axes.grid(onoff, "minor", color="#e0e0e0", linestyle=":")
            axes.tick_params(which="both", direction="in")
        else:
            axes.grid(onoff, "major")
            axes.grid(onoff, "minor")
    else:
        pp.minorticks_on()
        if onoff:
            pp.grid(onoff, "major", color="#d8d8d8")
            pp.grid(onoff, "minor", color="#e0e0e0", linestyle=":")
            pp.tick_params(which="both", direction="in")
        else:
            pp.grid(onoff, "major")
            pp.grid(onoff, "minor")

def set_color_cycle(c):
    if c in color_cycles.keys():
        pp.rcParams["axes.prop_cycle"] = cycler(color=color_cycles[c])
    else:
        print(f"plotdefs.set_color_cycle: cycle name '{c}' undefined. Colors were not modified.")

def opacity(i, N, m=0.2, p=2):
    return m + (1 - m) * (i / N)**p

def dot_dash(n):
    return n * [1, 1] + [8, 1]

class Slicer:
    def __init__(self):
        return

    def __getitem__(self, slice):
        return slice

S = Slicer()

class Plots:
    def __init__(self, plots):
        assert all(p.fig is plots[0].fig for p in plots)
        self.fig = plots[0].fig
        self.plots = plots

    def __getitem__(self, pos):
        return self.plots[pos]

    def tight_layout(self, *args, **kwargs):
        X = self.fig.tight_layout(*args, **kwargs)
        self.outputs.append(X)
        return self

    def savefig(self, *args, **kwargs):
        X = self.fig.savefig(*args, **kwargs)
        self.outputs.append(X)
        return self

    def show(self):
        pp.show()
        return self

    def close(self):
        pp.close(self.fig)

    def f(self, f, *args, **kwargs):
        X = f(*args, **kwargs)
        self.outputs.append(X)
        return self

class Plotter:
    def __init__(self, fig=None, ax=None):
        if fig is None or ax is None:
            self.fig, self.ax = pp.subplots()
        else:
            self.fig = fig
            self.ax = ax
        self.outputs = list()
        self.im = None
        self.cbar = None

    @staticmethod
    def new(*args, **kwargs):
        if kwargs.get("fig", None) is None \
                or kwargs.get("ax", None) is None:
            fig, ax = pp.subplots(*args, **kwargs)
        else:
            fig, ax = kwargs["fig"], kwargs["ax"]
        if isinstance(ax, (np.ndarray, list, tuple)):
            return [Plotter(fig=fig, ax=a) for a in ax]
        else:
            return Plotter(fig=fig, ax=ax)

    @staticmethod
    def new_3d(*args, **kwargs):
        fig = pp.figure(*args, **kwargs)
        ax = p3.Axes3D(fig)
        return Plotter(fig=fig, ax=ax)

    @staticmethod
    def new_gridspec(gridspec_kw, pos, *args, **kwargs):
        fig = kwargs.get("fig", pp.figure(*args, **kwargs))
        gs = fig.add_gridspec(**gridspec_kw)
        ax = [fig.add_subplot(gs[p]) for p in plots]
        return [Plotter(fig=fig, ax=a) for a in ax]

    def twinx(self, *args, **kwargs):
        return (
            self,
            Plotter(fig=self.fig, ax=self.ax.twinx())
        )

    def twiny(self, *args, **kwargs):
        return (
            self,
            Plotter(self.fig, ax=self.ax.twiny())
        )

    def sharex(self, *args, **kwargs):
        X = self.ax.sharex(*args, **kwargs)
        self.outputs.append(X)
        return self

    def sharey(self, *args, **kwargs):
        X = self.ax.sharey(*args, **kwargs)
        self.outputs.append(X)
        return self

    def plot(self, *args, **kwargs):
        X = self.ax.plot(*args, **kwargs)
        self.outputs.append(X)
        return self

    def plot_surface(self, *args, **kwargs):
        X = self.ax.plot_surface(*args, **kwargs)
        self.outputs.append(X)
        return self

    def plot_trisurf(self, *args, **kwargs):
        X = self.ax.plot_trisurf(*args, **kwargs)
        self.outputs.append(X)
        return self

    def errorbar(self, *args, **kwargs):
        X = self.ax.errorbar(*args, **kwargs)
        self.outputs.append(X)
        return self

    def semilogx(self, *args, **kwargs):
        X = self.ax.semilogx(*args, **kwargs)
        self.outputs.append(X)
        return self

    def semilogy(self, *args, **kwargs):
        X = self.ax.semilogy(*args, **kwargs)
        self.outputs.append(X)
        return self

    def loglog(self, *args, **kwargs):
        X = self.ax.loglog(*args, **kwargs)
        self.outputs.append(X)
        return self

    def scatter(self, *args, **kwargs):
        X = self.ax.scatter(*args, **kwargs)
        self.outputs.append(X)
        return self

    def contour(self, *args, **kwargs):
        mut = kwargs.pop("mut") if "mut" in kwargs.keys() else True
        X = self.ax.contour(*args, **kwargs)
        self.outputs.append(X)
        if mut:
            self.im = X
        return self

    def contourf(self, *args, **kwargs):
        mut = kwargs.pop("mut") if "mut" in kwargs.keys() else True
        X = self.ax.contourf(*args, **kwargs)
        self.outputs.append(X)
        if mut:
            self.im = X
        return self

    def axhline(self, *args, **kwargs):
        X = self.ax.axhline(*args, **kwargs)
        self.outputs.append(X)
        return self

    def axvline(self, *args, **kwargs):
        X = self.ax.axvline(*args, **kwargs)
        self.outputs.append(X)
        return self

    def axline(self, *args, **kwargs):
        X = self.ax.axline(*args, **kwargs)
        self.outputs.append(X)
        return self

    def fill(self, *args, **kwargs):
        X = self.ax.fill(*args, **kwargs)
        self.outputs.append(X)
        return self

    def fill_between(self, *args, **kwargs):
        X = self.ax.fill_between(*args, **kwargs)
        self.outputs.append(X)
        return self

    def imshow(self, *args, **kwargs):
        mut = kwargs.pop("mut") if "mut" in kwargs.keys() else True
        X = self.ax.imshow(*args, **kwargs)
        self.outputs.append(X)
        if mut:
            self.im = X
        return self

    def hist(self, *args, **kwargs):
        X = self.ax.hist(*args, **kwargs)
        self.outputs.append(X)
        return self

    def hist2d(self, *args, **kwargs):
        X = self.ax.hist2d(*args, **kwargs)
        self.outputs.append(X)
        return self

    def bar(self, *args, **kwargs):
        X = self.ax.bar(*args, **kwargs)
        self.outputs.append(X)
        return self

    def quiver(self, *args, **kwargs):
        X = self.ax.quiver(*args, **kwargs)
        self.outputs.append(X)
        return self

    def streamplot(self, *args, **kwargs):
        X = self.ax.streamplot(*args, **kwargs)
        self.outputs.append(X)
        return self

    def arrow(self, *args, **kwargs):
        X = self.ax.arrow(*args, **kwargs)
        self.outputs.append(X)
        return self

    def indicate_inset(self, *args, **kwargs):
        X = self.ax.indicate_inset(*args, **kwargs)
        self.outputs.append(X)
        return self

    def indicate_inset_zoom(self, *args, **kwargs):
        X = self.ax.indicate_inset_zoom(*args, **kwargs)
        self.outputs.append(X)
        return self

    def inset_axes(self, *args, **kwargs):
        mut = kwargs.pop("mut") if "mut" in kwargs.keys() else True
        X = self.ax.inset_axes(*args, **kwargs)
        self.outputs.append(X)
        if mut:
            self.ax = X
        return self

    def secondary_xaxis(self, *args, **kwargs):
        mut = kwargs.pop("mut") if "mut" in kwargs.keys() else True
        X = self.ax.secondary_xaxis(*args, **kwargs)
        self.outputs.append(X)
        if mut:
            self.ax = X
        return self

    def secondary_yaxis(self, *args, **kwargs):
        mut = kwargs.pop("mut") if "mut" in kwargs.keys() else True
        X = self.ax.secondary_yaxis(*args, **kwargs)
        self.outputs.append(X)
        if mut:
            self.ax = X
        return self

    def text(self, *args, **kwargs):
        X = self.ax.text(*args, **kwargs)
        self.outputs.append(X)
        return self

    def get_xlim(self, *args, **kwargs):
        return self.ax.get_xlim(*args, **kwargs)

    def get_ylim(self, *args, **kwargs):
        return self.ax.get_ylim(*args, **kwargs)

    def get_clim(self, *args, **kwargs):
        return self.im.get_clim(*args, **kwargs)

    def get_xticks(self, *args, **kwargs):
        return self.ax.get_xticks(*args, **kwargs)

    def get_yticks(self, *args, **kwargs):
        return self.ax.get_yticks(*args, **kwargs)

    def get_cticks(self, *args, **kwargs):
        return self.cbar.get_ticks(*args, **kwargs)

    def get_xticklabels(self, *args, **kwargs):
        return self.ax.get_xticklabels(*args, **kwargs)

    def get_yticklabels(self, *args, **kwargs):
        return self.ax.get_yticklabels(*args, **kwargs)

    def get_cticklabels(self, *args, **kwargs):
        return self.cbar.get_ticklabels(*args, **kwargs)

    def get_xlabel(self, *args, **kwargs):
        return self.ax.get_xlabel(*args, **kwargs)

    def get_ylabel(self, *args, **kwargs):
        return self.ax.get_ylabel(*args, **kwargs)

    def get_clabel(self, *args, **kwargs):
        return self.cbar.get_label(*args, **kwargs)

    def get_title(self, *args, **kwargs):
        return self.ax.get_title(*args, **kwargs)

    def set_xlim(self, *args, **kwargs):
        X = self.ax.set_xlim(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_ylim(self, *args, **kwargs):
        X = self.ax.set_ylim(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_clim(self, *args, **kwargs):
        X = self.im.set_clim(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_xticks(self, *args, **kwargs):
        X = self.ax.set_xticks(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_yticks(self, *args, **kwargs):
        X = self.ax.set_yticks(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_cticks(self, *args, **kwargs):
        X = self.cbar.set_ticks(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_xticklabels(self, *args, **kwargs):
        X = self.ax.set_xticklabels(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_yticklabels(self, *args, **kwargs):
        X = self.ax.set_yticklabels(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_cticklabels(self, *args, **kwargs):
        X = self.cbar.set_cticklabels(*args, **kwargs)
        self.outputs.append(X)
        return self

    def tick_params(self, *args, **kwargs):
        X = self.ax.tick_params(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_xlabel(self, *args, **kwargs):
        X = self.ax.set_xlabel(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_ylabel(self, *args, **kwargs):
        X = self.ax.set_ylabel(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_clabel(self, *args, **kwargs):
        X = self.cbar.set_label(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_title(self, *args, **kwargs):
        X = self.ax.set_title(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set(self, **kwargs):
        which = kwargs.get("which", "ax")
        if which == "ax":
            O = self.ax
        elif which == "fig":
            O = self.fig
        elif which == "cbar":
            O = self.cbar
        elif which == "im":
            O = self.im
        else:
            raise Exception("invalid 'which'")
        X = O.set(**kwargs)
        self.outputs.append(X)
        return self

    def invert_xaxis(self, *args, **kwargs):
        X = self.ax.invert_xaxis(*args, **kwargs)
        self.outputs.append(X)
        return self

    def invert_yaxis(self, *args, **kwargs):
        X = self.ax.invert_yaxis(*args, **kwargs)
        self.outputs.append(X)
        return self

    def colorbar(self, *args, **kwargs):
        mut = kwargs.pop("mut") if "mut" in kwargs.keys() else True
        X = self.fig.colorbar(self.im, *args, **kwargs)
        self.outputs.append(X)
        if mut:
            self.cbar = X
        return self

    def grid(self, *args, **kwargs):
        X = self.ax.grid(*args, **kwargs)
        self.outputs.append(X)
        return self

    def ggrid(self, onoff=True, *args, **kwargs):
        X = grid(onoff, self.ax)
        self.outputs.append(X)
        return self

    def legend(self, *args, **kwargs):
        X = self.ax.legend(*args, **kwargs)
        self.outputs.append(X)
        return self

    def tight_layout(self, *args, **kwargs):
        X = self.fig.tight_layout(*args, **kwargs)
        self.outputs.append(X)
        return self

    def set_box_aspect(self, *args, **kwargs):
        X = self.ax.set_box_aspect(*args, **kwargs)
        self.outputs.append(X)
        return self

    def axis(self, *args, **kwargs):
        X = self.ax.axis(*args, **kwargs)
        self.outputs.append(X)
        return self

    def savefig(self, *args, **kwargs):
        X = self.fig.savefig(*args, **kwargs)
        self.outputs.append(X)
        return self

    def show(self):
        pp.show()
        return self

    def close(self):
        pp.close(self.fig)

    def f(self, f, *args, **kwargs):
        X = f(*args, **kwargs)
        self.outputs.append(X)
        return self

class FigSize:
    def __init__(self, wh):
        assert isinstance(wh, (list, tuple))
        assert len(wh) == 2
        self.wh = list(wh)

    def _opcheck(self, other):
        assert isinstance(other, (int, float, list, tuple, FigSize))
        if isinstance(other, (list, tuple)):
            assert len(other) == 2

    def __abs__(self, /):
        return FigSize([abs(self.__w), abs(self.__h)])

    def __pos__(self, /):
        return self
    
    def __neg__(self, /):
        return FigSize([-self.__w, -self.__h])

    def __invert__(self, /):
        return FigSize([self.__h, self.__w])

    def __contains__(self, val, /):
        return val in self.__wh

    def __eq__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return (self.__w == val) and (self.__h == val)
        elif isinstance(val, (list, tuple)):
            return (self.__w == val[0]) and (self.__h == val[0])
        elif isinstance(val, FigSize):
            return (self.__w == val.w) and (self.__h == val.h)

    def __ne__(self, val, /):
        return not (self == val)

    def __getitem__(self, key, /):
        assert key in [0, "w", 1, "h"]
        if key in [0, "w"]:
            return self.__w
        elif key in [1, "h"]:
            return self.__h

    def __setitem__(self, key, val, /):
        assert key in [0, "w", 1, "h"]
        if key in [0, "w"]:
            self.w = val
        elif key in [1, "h"]:
            self.h = val

    def __iter__(self, /):
        return iter(self.__wh)

    def __reversed__(self, /):
        return reversed(self.__wh)

    def __len__(self, /):
        return len(self.__wh)

    @property # wh
    def wh(self):
        return self.__wh
    @wh.setter
    def wh(self, wh, /):
        assert isinstance(wh, (list, tuple))
        assert len(wh) == 2
        assert isinstance(wh[0], (int, float)) and isinstance(wh[1], (int, float))
        self.__wh = list(wh)
        self.__w = wh[0]
        self.__h = wh[1]

    @property # w
    def w(self):
        return self.__w
    @w.setter
    def w(self, w, /):
        assert isinstance(w, (int, float))
        self.__wh[0] = w
        self.__w = w

    @property # h
    def h(self):
        return self.__h
    @h.setter
    def h(self, h, /):
        assert isinstance(h, (int, float))
        self.__wh[1] = h
        self.__h = h

    def __add__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([self.__w+val, self.__h+val])
        elif isinstance(val, (list, tuple)):
            return FigSize([self.__w+val[0], self.__h+val[1]])
        elif isinstance(val, FigSize):
            return FigSize([self.__w+val.w, self.__h+val.h])
    def __radd__(self, val, /):
        return self.__add__(val)
    def __iadd__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            self.w = self.w + val
            self.h = self.h + val
        elif isinstance(val, (list, tuple)):
            self.w = self.w + val[0]
            self.h = self.h + val[1]
        elif isinstance(val, FigSize):
            self.w = self.w + val.w
            self.h = self.h + val.h

    def __sub__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([self.__w-val, self.__h-val])
        elif isinstance(val, (list, tuple)):
            return FigSize([self.__w-val[0], self.__h-val[1]])
        elif isinstance(val, FigSize):
            return FigSize([self.__w-val.w, self.__h-val.h])
    def __rsub__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([val-self.__w, val-self.__h])
        elif isinstance(val, (list, tuple)):
            return FigSize([val[0]-self.__w, val[1]-self.__h])
        elif isinstance(val, FigSize):
            return FigSize([val.w-self.__w, val.h-self.__h])
    def __isub__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            self.w = self.__w - val
            self.h = self.__h - val
        elif isinstance(val, (list, tuple)):
            self.w = self.__w - val[0]
            self.h = self.__h - val[1]
        elif isinstance(val, FigSize):
            self.w = self.__w - val.w
            self.h = self.__h - val.h

    def __mul__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([self.__w*val, self.__h*val])
        elif isinstance(val, (list, tuple)):
            return FigSize([self.__w*val[0], self.__h*val[1]])
        elif isinstance(val, FigSize):
            return FigSize([self.__w*val.w, self.__h*val.h])
    def __rmul__(self, val, /):
        return self.__mul__(val)
    def __isub__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            self.w = self.__w * val
            self.h = self.__h * val
        elif isinstance(val, (list, tuple)):
            self.w = self.__w * val[0]
            self.h = self.__h * val[1]
        elif isinstance(val, FigSize):
            self.w = self.__w * val.w
            self.h = self.__h * val.h

    def __truediv__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([self.__w/val, self.__h/val])
        elif isinstance(val, (list, tuple)):
            return FigSize([self.__w/val[0], self.__h/val[1]])
        elif isinstance(val, FigSize):
            return FigSize([self.__w/val.w, self.__h/val.h])
    def __rtruediv__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([val/self.__w, val/self.__h])
        elif isinstance(val, (list, tuple)):
            return FigSize([val[0]/self.__w, val[1]/self.__h])
        elif isinstance(val, FigSize):
            return FigSize([val.w/self.__w, val.h/self.__h])
    def __itruediv__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            self.w = self.__w / val
            self.h = self.__h / val
        elif isinstance(val, (list, tuple)):
            self.w = self.__w / val[0]
            self.h = self.__h / val[1]
        elif isinstance(val, FigSize):
            self.w = self.__w / val.w
            self.h = self.__h / val.h

    def __floordiv__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([self.__w//val, self.__h//val])
        elif isinstance(val, (list, tuple)):
            return FigSize([self.__w//val[0], self.__h//val[1]])
        elif isinstance(val, FigSize):
            return FigSize([self.__w//val.w, self.__h//val.h])
    def __rfloordiv__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([val//self.__w, val//self.__h])
        elif isinstance(val, (list, tuple)):
            return FigSize([val[0]//self.__w, val[1]//self.__h])
        elif isinstance(val, FigSize):
            return FigSize([val.w//self.__w, val.h//self.__h])
    def __ifloordiv__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            self.w = self.__w // val
            self.h = self.__h // val
        elif isinstance(val, (list, tuple)):
            self.w = self.__w // val[0]
            self.h = self.__h // val[1]
        elif isinstance(val, FigSize):
            self.w = self.__w // val.w
            self.h = self.__h // val.h

    def __mod__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([self.__w%val, self.__h%val])
        elif isinstance(val, (list, tuple)):
            return FigSize([self.__w%val[0], self.__h%val[1]])
        elif isinstance(val, FigSize):
            return FigSize([self.__w%val.w, self.__h%val.h])
    def __rmod__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([val%self.__w, val%self.__h])
        elif isinstance(val, (list, tuple)):
            return FigSize([val[0]%self.__w, val[1]%self.__h])
        elif isinstance(val, FigSize):
            return FigSize([val.w%self.__w, val.h%self.__h])
    def __imod__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            self.w = self.__w % val
            self.h = self.__h % val
        elif isinstance(val, (list, tuple)):
            self.w = self.__w % val[0]
            self.h = self.__h % val[1]
        elif isinstance(val, FigSize):
            self.w = self.__w % val.w
            self.h = self.__h % val.h

    def __pow__(self, val, mod=None, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            return FigSize([pow(self.__w, val, mod), pow(self.__h, val, mod)])
        elif isinstance(val, (list, tuple)):
            return FigSize([pow(self.__w, val[0], mod), pow(self.__h, val[1], mod)])
        elif isinstance(val, FigSize):
            return FigSize([pow(self.__w, val.w, mod), pow(self.__h, val.h, mod)])
    def __ipow__(self, val, /):
        self._opcheck(val)
        if isinstance(val, (int, float)):
            self.w = self.__w ** val
            self.h = self.__h ** val
        elif isinstance(val, (list, tuple)):
            self.w = self.__w ** val[0]
            self.h = self.__h ** val[1]
        elif isinstance(val, FigSize):
            self.w = self.__w ** val.w
            self.h = self.__h ** val.h

    def __repr__(self, /):
        return "FigSize("+str(self.__wh)+")"

    def __str__(self, /):
        return "FigSize("+str(self.__wh)+")"

