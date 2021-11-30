import networkx as nx
import numpy as np
import PySimpleGUI as sg
import time

def _scale_coords(coords, canvas_center, padding=0.12, move=(0., 0.)):
    """
    Rescale 2D coordinates based on canvas center and padding.
    """
    out_coords = list()
    for coord, center, offset in zip(coords, canvas_center, move):
        scale_factor = center - padding * center
        out_coords.append(center+coord*scale_factor-offset)
    return tuple(out_coords)

def draw_graph(canvas,
               graph,
               canvas_center,
               gen_coords=True,
               move=(0,0),
               padding=0.12,
               radius_scale=80,
               colors={},
               methods={},
               coordinates=None):
    """
    Draw a graph anew on an existing canvas object.
    To zoom in add a negative padding to the coordinates and to zoom out
    a higher than default padding.
    """
    # erease old canvas
    canvas.erase()
    # get the layout
    if gen_coords:
        coord_dict = nx.kamada_kawai_layout(graph, pos=coordinates)
    else:
        coord_dict = coordinates
    radius = 1/np.sqrt(len(graph.nodes)) * radius_scale
    # draw lines connecting nodes
    for ndxA, ndxB in graph.edges():
        coordA = _scale_coords(coord_dict[ndxA], canvas_center, padding=padding, move=move)
        coordB = _scale_coords(coord_dict[ndxB], canvas_center, padding=padding, move=move)
        canvas.draw_line(point_from=coordA, point_to=coordB, color='black')
    # draw node pictograms
    for node, coord in coord_dict.items():
        location = _scale_coords(coord, canvas_center, padding=padding, move=move)

        # use custom node color
        try:
            color=colors[graph.nodes[node]["resname"]]
        except KeyError:
            color='gray'

        # use custom node drawing method
        try:
           method = getattr(canvas, methods[graph.nodes[node]["resname"]])
           method(radius=radius, center_location=location, fill_color=color)
        except KeyError:
           canvas.draw_circle(radius=radius, center_location=location, fill_color=color)
    return coord_dict
