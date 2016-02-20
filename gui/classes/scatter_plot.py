from enthought.traits.api import HasTraits, Instance, Int
from enthought.traits.ui.api import View, Group, Item
from enthought.enable.api import ColorTrait
from enthought.enable.component_editor import ComponentEditor
from enthought.chaco.api import marker_trait, Plot, ArrayPlotData

class ScatterPlot(HasTraits):

    plot = Instance(Plot)
    color = ColorTrait("blue")
    marker = marker_trait
    marker_size = Int(4)

    traits_view1 = View(
            Group (
                Group(Item('color', label="Color", style="custom"),
                    Item('marker', label="Marker"),
                    Item('marker_size', label="Size"),
                    orientation = "vertical"),
                Group(Item('plot', editor=ComponentEditor(), show_label=False)),
                orientation = "horizontal"),
            width=800, height=600, resizable=True, title="Chaco Plot")

    def __init__(self, x, y, color = "blue"):
        super(ScatterPlot, self).__init__()
        plotdata = ArrayPlotData(x = x, y = y)
        plot = Plot(plotdata)
        self.renderer = plot.plot(("x", "y"), type = "scatter", color = color)[0]
        self.plot = plot
        self.color = color
        self.configure_traits()

    def _color_changed(self):
        self.renderer.color = self.color

    def _marker_changed(self):
        self.renderer.marker = self.marker

    def _marker_size_changed(self):
        self.renderer.marker_size = self.marker_size
