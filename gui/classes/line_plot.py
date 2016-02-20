from enthought.traits.api import HasTraits, Instance, Int
from enthought.traits.ui.api import View, Group, Item
from enthought.enable.api import ColorTrait
from enthought.enable.component_editor import ComponentEditor
from enthought.chaco.api import Plot, ArrayPlotData

class LinePlot(HasTraits):

    plot = Instance(Plot)
    color = ColorTrait("blue")
    line_width = Int(4)

    traits_view1 = View(
            Group (
                Group(Item('color', label="Color", style="custom"),
                    Item('line_width', label="Width"),
                    orientation = "vertical"),
                Group(Item('plot', editor=ComponentEditor(), show_label=False)),
                orientation = "horizontal"),
            width=800, height=600, resizable=True, title="Chaco Plot")

    def __init__(self, x, y, color = "blue", line_width = 4):
        super(LinePlot, self).__init__()
        plotdata = ArrayPlotData(x = x, y = y)
        plot = Plot(plotdata)
        self.renderer = plot.plot(("x", "y"), type = "line", color = color, line_width = line_width)[0]
        self.plot = plot
        self.color = color
        self.line_width = line_width
        self.configure_traits()

    def _color_changed(self):
        self.renderer.color = self.color

    def _line_width_changed(self):
        self.renderer.line_width = self.line_width
