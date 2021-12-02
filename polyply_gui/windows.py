import PySimpleGUI as sg
import os.path
import json
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph
import polyply
from polyply import DATA_PATH
from polyply.src.load_library import load_library
from polyply.src.gen_seq import _branched_graph
from . import IMG_PATH

class WindowCreator():
    """
    Base class for creating sg windows.
    """
    def __init__(self, modal=True, finalize=True, **kwargs):
        self.modal = modal
        self.finalize = finalize
        self.kwargs = kwargs

    def create_window(self):
        """
        Creates new instance of :class:`sg.Window`.
        """
        new_window = sg.Window(modal=self.modal,
                               finalize=self.finalize, **self.kwargs)
        return new_window


class ChainArchitechtureWindow(WindowCreator):
    """
    Windows handling the I/O for chain structure
    of blocks.
    """

    def __init__(self,
                 combo1_title="type",
                 combo1_values=list(),
                 combo1_visible=True,
                 in1_title="repeat units",
                 in1_value="1",
                 in1_visible=True,
                 in2_title="degree of branching",
                 in2_value="1",
                 in2_visible=False,
                 tacticity=False,
                 **kwargs):
        """
        Initalize the object with variable input for the columns.
        """
        rows = [[sg.Text("{:<15s}".format(combo1_title), visible=combo1_visible)],
                [sg.Combo(combo1_values,
                          size=(40, 4),
                          enable_events=False,
                          key='-MONOMERS-',
                          visible=combo1_visible)],
                 [sg.Text(in1_title, visible=in1_visible)],
                 [sg.In(enable_events=False,
                             key='-NMON-',
                          size=(40, 4),
                          visible=in1_visible,
                          default_text=in1_value)],
                  [sg.Text(in2_title, visible=in2_visible)],
                  [sg.In(enable_events=False,
                         key='-BRANCH-',
                         size=(40, 4),
                         visible=in2_visible)],
                  [sg.Text("set tacticity", visible=tacticity)],
                  [sg.Listbox(values=['atactic', 'isotactic-R', 'isotactic-S', 'syndiotactic'],
                              enable_events=False, key='tacticity', visible=tacticity)],
                  [sg.Button('add', enable_events=True, key='add_block')]]
        layout = [[sg.Column(rows)]]
        super().__init__(layout=layout, **kwargs)

class GraphViewerWindow(WindowCreator):
    """
    Window for the graph viewer element. This class also implements
    zoom related methods.
    """

    def __init__(self, canvas_size, background_color='white', **kwargs):
        self.canvas_size = canvas_size
        self.background_color = background_color
        rows = [[sg.Graph(background_color=background_color,
                          canvas_size=canvas_size,
                          graph_bottom_left=(0,0),
                          graph_top_right=canvas_size,
                          key='graph_event',
                          drag_submits=True,
                          enable_events=True)],
                [sg.Text("Zoom"), sg.Button('+', enable_events=True, key='zoom_in'), sg.Button('-', enable_events=True, key='zoom_out')]]

        layout = [[sg.Column(rows)]]
        super().__init__(layout=layout, **kwargs)

class AddConnectionWindow(WindowCreator):
    """
    Windows handling the I/O for chain structure
    of blocks.
    """

    def __init__(self,
                 combo1_title="type",
                 combo1_values=list(),
                 combo1_visible=True,
                 in1_title="repeat units",
                 in1_key="connect_blocks",
                 in1_visible=True,
                 tacticity=False,
                 **kwargs):
        """
        Initalize the object with variable input for the columns.
        """
        print(combo1_visible)
        rows = [[sg.Text("block"), sg.In(key='idA', size=(5, 1)), sg.Text("resid"), sg.In(key="residA", size=(5, 1))],
                [sg.Text("block"), sg.In(key='idB', size=(5, 1)), sg.Text("resid"), sg.In(key="residB", size=(5, 1))],
                [sg.Button(in1_title, enable_events=True, key=in1_key)],
                [sg.Text("select linkage type", visible=combo1_visible)]]
               #[sg.Listbox(combo1_values,
               #          size=(40, 4),
               #          enable_events=False,
               #          key='edge_label',
               #          visible=combo1_visible)]]

        layout = [[sg.Column(rows)]]
        super().__init__(layout=layout, **kwargs)

class MainWindow(WindowCreator):
    """
    Polyply GUI main window.
    """
    def __init__(self, libs, **kwargs):
        add_blocks_rows = [[sg.Text("Add blocks", justification="center", size=(20, 1), font='bold')],
                           [sg.Text("{:<15s}".format("force-field")),],
                           [sg.Combo(libs, size=(20, 2), enable_events=True, key='set_force_field')],
                           [sg.Text("{:<15s}".format("chain architecture")),],
                           [sg.Button(image_filename=IMG_PATH+"/linear_button.png", size=(20, 2), image_size=(200, 120),
                                      image_subsample=4, enable_events=True, key='linear_block', button_color=(None, "white"),
                                      mouseover_colors=("blue", "blue")),
                            sg.Button(image_filename=IMG_PATH+"/tree_button.png", size=(20, 2), image_size=(200, 120),
                                      image_subsample=6, enable_events=True, key='tree_block', button_color=(None, "white"),
                                      mouseover_colors=("blue", "blue"))],
                            [sg.Text("{:<26s}".format("linear")), sg.Text("tree")],
                            [sg.Button(image_filename=IMG_PATH+"/comb.png", size=(20, 2), image_size=(200, 120),
                                       image_subsample=4, enable_events=True, key='comb_block', button_color=(None, "white"),
                                       mouseover_colors=("blue", "blue")),
                             sg.Button(image_filename=IMG_PATH+"/network.png", size=(20, 2), image_size=(200, 120),
                                       image_subsample=4, enable_events=True, key='carbohydrate_block', button_color=(None, "white"),
                                       mouseover_colors=("blue", "blue"))],
                            [sg.Text("{:<26s}".format("comb")), sg.Text("carbohydrate")],
                            [sg.Text("From File")],
                            [sg.In(size=(25, 1), enable_events=True, key="load_file"),sg.FileBrowse()]]

        connect_column = [[sg.Text("Modify blocks", justification="center", size=(20, 1), font='bold')],
                          [sg.Text("Current blocks")],
                          [sg.Listbox(values=[], size=(20, 5), key="-SEQ-")],
                          [sg.Button('delete selected block', enable_events=True, key='remove_block')],
                          [sg.Button(image_filename=IMG_PATH+"/staple.png", size=(10, 5), image_size=(150, 200),
                                      image_subsample=2, enable_events=True, key='open_connect_window', button_color=(None, "white"),
                                      mouseover_colors=("blue", "blue")),
                           sg.Button(image_filename=IMG_PATH+"/no_staple.png", size=(10, 5), image_size=(150, 200),
                                      image_subsample=3, enable_events=True, key='open_delete_window', button_color=(None, "white"),
                                      mouseover_colors=("blue", "blue")),
                            ],
                          [sg.Text("connect", size=(10, 1)), sg.Text("disconnect", size=(10, 1))],
                          #[sg.Text("block"), sg.In(key='idA', size=(5,1 )), sg.Text("resid"), sg.In(key="residA", size=(5, 1))],
                          #[sg.Text("block"), sg.In(key='idB', size=(5,1 )), sg.Text("resid"), sg.In(key="residB", size=(5, 1))],
                          #[sg.Button('add', enable_events=True, key='conncet_blocks'), sg.Button('delete', enable_events=True, key='remove_edge')]
                          ]

        sequence_viewer_column = [[sg.Text("I/O", justification="center", size=(20, 1), font='bold')],
                                  [sg.Text("First save the graph.", size=(20, 1))],
                                  [sg.Text("Then run generate itp.", size=(20, 1))],
                                  [sg.SaveAs(button_text="save graph", enable_events=True, key='write_seq_file')],
                                  [sg.SaveAs(button_text="generate itp file", enable_events=True, key='gen_itp')]]

        layout = [[sg.Column(add_blocks_rows, vertical_alignment='top'),
                   sg.VSeperator(),
                   sg.Column(connect_column, vertical_alignment='top'),
                   sg.VSeperator(),
                   sg.Column(sequence_viewer_column, vertical_alignment='top'),
                   ],
                   [sg.Text("Command Log")],
                   [sg.Listbox(enable_events=False,
                             key='update_log',
                             size=(85, 4),
                             visible=True,
                             values=["This is a polyply command"]
                             )]]

        super().__init__(layout=layout, **kwargs)
