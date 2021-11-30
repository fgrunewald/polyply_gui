import os
import PySimpleGUI as sg
import os.path
from pathlib import Path
import json
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph
import vermouth
from vermouth.gmx.itp_read import read_itp
from vermouth.gmx.gro import read_gro
from vermouth.pdb.pdb import read_pdb
from vermouth.graph_utils import make_residue_graph
import polyply
from polyply import DATA_PATH
from polyply.src.load_library import load_library
from polyply.src.gen_seq import _branched_graph, _random_replace_nodes_attribute
from .windows import ChainArchitechtureWindow, AddConnectionWindow
from .graph_drawing import draw_graph
from .polyply_runner import run_gen_itp

PARSERS = {"itp": read_itp,}

class EventHandler():
    """
    Class for handling the bookkeeping behind the different
    events in the polyply GUI. Each event is implemented as
    a class method, which take the window and the event
    that triggered them as input.
    """

    def __init__(self, base_window, graph_viewer, canvas_size):
        self.graph = nx.Graph()
        self.blocks = []
        self.base_window = base_window
        self.graph_viewer = graph_viewer
        self.force_field = None
        self.monomers = []
        self.seq_list = []
        self.canvas_size = canvas_size
        self.canvas_center = (canvas_size[0]/2., canvas_size[1]/2.)
        self.zoom_factor = 0
        self.seq_path = None
        self.itp_path = None
        self.coords = None
        self.arch_args = {"tree_block": {'title':'tree block' ,
                                     'combo1_title': "monomer type",
                                     'combo1_values': list(),
                                     'combo1_visible': True,
                                     'in1_title': "generations",
                                     'in1_value': "1",
                                     'in1_visible': True,
                                     'in2_title': "degree of branching",
                                     'in2_value': "1",
                                     'in2_visible': True,},
                          "comb_block": {'title':'comb block',
                                     'combo1_title': "monomer type",
                                     'combo1_values': list(),
                                     'combo1_visible': True,
                                     'in1_title': "repeat units",
                                     'in1_value': "1",
                                     'in1_visible': True,
                                     'in2_title': "nnn",
                                     'in2_value': "1",
                                     'in2_visible': False},
                          "carbohydrate_block":{'title': 'monomer',
                                   'combo1_title': "monomer type",
                                   'combo1_values': list(),
                                   'combo1_visible': True,
                                   'in1_title': "",
                                   'in1_value': "",
                                   'in1_visible': False,
                                   'in2_title': "",
                                   'in2_value': "",
                                   'in2_visible': False},
                          "linear_block": {'title': 'linear block',
                                   'combo1_title': "monomer type",
                                   'combo1_values': list(),
                                   'combo1_visible': True,
                                   'in1_title': "repeat units",
                                   'in1_value': "1",
                                   'in1_visible': True,
                                   'in2_title': "",
                                   'in2_value': "",
                                   'in2_visible': False,}}

    def set_force_field(self, window, event, values):
        self.force_field = window["set_force_field"].get()
        library = load_library("libs", [self.force_field], [])
        self.monomers = list(library.blocks.keys())

    def open_architecture_window(self, window, event, values):
        window_args = self.arch_args[event]
        if self.force_field not in ['martini3', 'martini2']:
            window_args['tacticity'] = True

        new_window = ChainArchitechtureWindow(**self.arch_args[event]).create_window()
        new_window["-MONOMERS-"].update(values=self.monomers)

        if event != "tree_block":
            new_window["-BRANCH-"].update(value=1)

        if event == "carbohydrate_block":
            new_window["-NMON-"].update(value=1)

        return new_window

    def linear_block(self, window, event, values):
        self.open_architecture_window(window, event, values)

    def tree_block(self, window, event, values):
        self.open_architecture_window(window, event, values)

    def comb_block(self, window, event, values):
        self.open_architecture_window(window, event, values)

    def carbohydrate_block(self, window, event, values):
        self.open_architecture_window(window, event, values)

    def open_connect_window(self, window, event, values):
        new_window = AddConnectionWindow(title="add links", in1_title="add", in1_key="conncet_blocks").create_window()

    def open_delete_window(self, window, event, values):
        new_window = AddConnectionWindow(title="remove links",
                                         combo1_visible=False,
                                         in1_title="remove",
                                         in1_key="remove_block").create_window()

    def add_block(self, window, event, values):
        block_count = len(self.blocks) + 1
        # get the user input
        n_mon = values['-NMON-']
        branching_f = values['-BRANCH-']
        monomer = window['-MONOMERS-'].get()
        # update the sequence list
        seq_element = " ".join([str(block_count), monomer, n_mon])
        self.seq_list.append(seq_element)
        self.base_window["-SEQ-"].update(self.seq_list)
        # update nx-graph
        last = len(self.graph.nodes)
        new_graph = _branched_graph(monomer, int(branching_f), int(n_mon))
        # set the tacticity attribute
        if self.force_field not in ['martini3', 'martini2']:
            tacticity = window['tacticity'].get()[0]
            if tacticity == 'atactic':
                print(tacticity)
                self.graph = _random_replace_nodes_attribute(self.graph,
                                                             attribute_values=["R", "S"],
                                                             weights=[0.5, 0.5],
                                                             attribute='tacticity',
                                                             nodes=[],
                                                             seed=None)
            elif tacticity == 'isotactic-R':
                nx.set_node_attributes(self.graph, 'R', 'tacticity')
            elif tacticity == 'isotactic-S':
                nx.set_node_attributes(self.graph, 'S', 'tacticity')

        self.graph = nx.disjoint_union(self.graph, new_graph)
        self.blocks.append(range(last, len(self.graph)))
        # draw the residue graph
        canvas = self.graph_viewer['graph_event']
        self.coords = draw_graph(canvas, self.graph, canvas_center=self.canvas_center)

    def conncet_blocks(self, window, event, values):
        idxA = int(values["idA"]) - 1
        idxB = int(values["idB"]) - 1
        residA = int(values["residA"]) - 1
        residB = int(values["residB"]) - 1
        nodeA = self.blocks[idxA][residA]
        nodeB = self.blocks[idxB][residB]
        self.graph.add_edge(nodeA, nodeB)
        canvas = self.graph_viewer['graph_event']
        self.coords = draw_graph(canvas, self.graph, canvas_center=self.canvas_center, coordinates=self.coords)

    def remove_edge(self, window, event, values):
        idxA = int(values["idA"])
        idxB = int(values["idB"])
        residA = int(values["residA"])
        residB = int(values["residB"])
        nodeA = blocks[idxA][residA]
        nodeB = blocks[idxB][residB]
        graph.remove_edge(nodeA, nodeB)
        canvas = self.graph_viewer['graph_event']
        self.coords = draw_graph(canvas, graph, canvas_center=self.canvas_center, coordinates=self.coords)

    def remove_block(self, window, event, values):
        seq_idx = int(window["-SEQ-"].get()[0].split()[0]) - 1
        nodes  = self.blocks[seq_idx]
        for node in nodes:
            self.graph.remove_node(node)
        del self.blocks[seq_idx]
        del self.seq_list[seq_idx]
        for idx, item in enumerate(self.seq_list):
            jdx, monomer, n_mon = item.split()
            self.seq_list[idx] = " ".join([str(idx), monomer, n_mon])
        window["-SEQ-"].update(self.seq_list)
        canvas = self.graph_viewer['graph_event']
        self.coords = draw_graph(canvas, self.graph, canvas_center=self.canvas_center, coordinates=self.coords)

    def write_seq_file(self, window, event, values):
        self.seq_path = values['write_seq_file']
        if self.seq_path:
            g_json = json_graph.node_link_data(self.graph)
            with open(self.seq_path, "w") as file_handle:
                json.dump(g_json, file_handle, indent=2)

    def zoom_in(self, window, event, values):
        canvas = self.graph_viewer['graph_event']
        self.zoom_factor += -0.12
        self.coords = draw_graph(canvas, self.graph, gen_coords=False,
                   canvas_center=self.canvas_center,
                   padding=self.zoom_factor, coordinates=self.coords)

    def zoom_out(self, window, event, values):
        canvas = self.graph_viewer['graph_event']
        self.zoom_factor += 0.12
        self.coords = draw_graph(canvas,
                                 self.graph,
                                 canvas_center=self.canvas_center,
                                 gen_coords=False,
                                 padding=self.zoom_factor,
                                 coordinates=self.coords)

    def graph_event(self, window, event, values):
        canvas = self.graph_viewer['graph_event']
        xloc, yloc = values['graph_event']
        deltax, deltay = self.canvas_center[0] - xloc, self.canvas_center[1] - yloc
        self.coords = draw_graph(canvas,
                                 self.graph,
                                 gen_coords=False,
                                 padding=self.zoom_factor,
                                 canvas_center=self.canvas_center,
                                 move=(deltax, deltay),
                                 coordinates=self.coords)

    def gen_itp(self, window, event, values):
        self.itp_path = values['gen_itp']
        if self.itp_path:
            output = run_gen_itp(self.seq_path, self.itp_path, self.force_field)
            values['update_log'] = output
            self.update_log(window, event, values)

    def update_log(self, window, event, values):
        log = window["update_log"].get()
        log += values["update_log"]
        window["update_log"].update(log)

    def load_file(self, window, event, values):
        path = Path(values['load_file'])
        if path.is_file():
            with open(path) as _file:
                lines = _file.readlines()
            file_extension = path.suffix.casefold()[1:]
            force_field = vermouth.forcefield.ForceField("dummy")
            PARSERS[file_extension](lines, force_field)
            mol_name = list(force_field.blocks.keys())[0]
            molecule = force_field.blocks[mol_name]
            molecule.make_edges_from_interaction_type('bonds')
            new_graph = make_residue_graph(molecule)
            last = len(self.graph.nodes)
            self.graph = nx.disjoint_union(self.graph, new_graph)
            self.blocks.append(range(last, len(self.graph)))
            # update sequence list
            block_count = len(self.blocks) + 1
            seq_element = " ".join([str(block_count), mol_name, str(len(new_graph.nodes))])
            self.seq_list.append(seq_element)
            self.base_window["-SEQ-"].update(self.seq_list)
            # draw the residue graph
            canvas = self.graph_viewer['graph_event']
            self.coords = draw_graph(canvas, self.graph, canvas_center=self.canvas_center, coordinates=self.coords)
