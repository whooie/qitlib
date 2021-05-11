"""
Handy matplotlib.pyplot settings.
"""
import matplotlib as mpl
import matplotlib.pyplot as pp
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.colors as cm
from cycler import cycler

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

color_cycles = {
    "matlab": [
        "#0072bd",
        "#d95319",
        "#edb120",
        "#7e2f8e",
        "#77ac30",
        "#4dbeee",
        "#a2142f"
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
        "#17becf"
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

bwr_signed_colors = [
    (0.000, "#101010"),
    (0.100, "#3f119d"),
    (0.350, "#3967d0"),
    (0.500, "#f0f0f0"),
    (0.625, "#f1b931"),
    (1.000, "#dd0000")
]
bwr_signed = cm.LinearSegmentedColormap.from_list("blue-white-red_signed", bwr_signed_colors)

bbw_unsigned_colors = [
    (0.000, "#000000"),
    (0.450, "#3b4568"),
    (0.600, "#586186"),
    (0.700, "#939cc4"),
    (1.000, "#ffffff")
]
bbw_unsigned = cm.LinearSegmentedColormap.from_list("black-blue-white_unsigned", bbw_unsigned_colors)

br_colors = [
    (0.000, "#101010"),
    (0.100, "#3967d0"),
    (1.000, "#dd0000")
]
br = cm.LinearSegmentedColormap.from_list("br", br_colors)

rcdefs = {
    "axes.grid"             : True,
    "axes.grid.which"       : "both",
    "axes.linewidth"        : 1.0,
    "axes.prop_cycle"       : cycler(color=color_cycles["best"]),
    "errorbar.capsize"      : 3.5,
    "figure.figsize"        : FigSize([14, 10.5]),
    "font.size"             : 24,
    "image.cmap"            : "bone",
    "legend.fontsize"       : "small",
    "lines.linewidth"       : 2.0,
    "lines.markersize"      : 8.0,
    "lines.markeredgewidth" : 3.5,
    "savefig.bbox"          : "tight",
    "savefig.pad_inches"    : 0.1,
    "text.latex.preamble"   : r"\usepackage{physics}\usepackage{siunitx}\usepackage{amsmath}",
    "xtick.direction"       : "in",
    "ytick.direction"       : "in",
}
for key in rcdefs:
    pp.rcParams[key] = rcdefs[key]

def figure3D(*fig_args, **fig_kwargs):
    fig = pp.figure(*fig_args, **fig_kwargs)
    ax = p3.Axes3D(fig)
    return fig, ax

def set_fonts_size(s):
    pp.rcParams["font.size"] = s

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
            axes.grid(onoff, "major", color="#b0b0b0")
            axes.grid(onoff, "minor", color="#b8b8b8", linestyle=":")
            axes.tick_params(which="both", direction="in")
        else:
            axes.grid(onoff, "major")
            axes.grid(onoff, "minor")
    else:
        pp.minorticks_on()
        if onoff:
            pp.grid(onoff, "major", color="#b0b0b0")
            pp.grid(onoff, "minor", color="#b8b8b8", linestyle=":")
            pp.tick_params(which="both", direction="in")
        else:
            pp.grid(onoff, "major")
            pp.grid(onoff, "minor")

def set_color_cycle(c):
    if c in color_cycles.keys():
        pp.rcParams["axes.prop_cycle"] = cycler(color=color_cycles[c])
    else:
        print(f"plotdefs.set_color_cycle: cycle name '{c}' undefined. Colors were not modified.")

def set_figsize(w, h):
    pp.rcParams["figure.figsize"] = FigSize([w, h])

def get_lims(ax=None):
    if ax:
        return ax.get_xlim(), ax.get_ylim()
    else:
        return pp.xlim(), pp.ylim()

def set_lims(xlim=None, ylim=None, ax=None):
    if ax:
        if xlim:
            ax.set_xlim(*xlim)
        if ylim:
            ax.set_ylim(*ylim)
    else:
        if xlim:
            pp.xlim(*xlim)
        if ylim:
            ax.ylim(*ylim)

def get_ticks(ax=None):
    if ax:
        return (ax.get_xticks(), ax.get_xticklabels()), (ax.get_yticks(), ax.get_yticklabels())
    else:
        return pp.xticks(), pp.yticks()

def set_ticks(xticks=None, yticks=None, ax=None):
    if ax:
        if xticks:
            if isinstance(xticks, tuple):
                ax.set_xticks(xticks[0])
                ax.set_xticklabels(xticks[1])
            else:
                ax.set_xticks(xticks)
        if yticks:
            if isinstance(yticks, tuple):
                ax.set_yticks(yticks[0])
                ax.set_yticklabels(yticks[1])
            else:
                ax.set_yticks(yticks)
    else:
        if xticks:
            if isinstance(xticks, tuple):
                pp.xticks(*xticks)
            else:
                pp.xticks(xticks)
        if yticks:
            if isinstance(yticks, tuple):
                pp.yticks(*yticks)
            else:
                pp.yticks(yticks)



