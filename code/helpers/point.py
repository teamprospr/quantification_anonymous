#!/usr/bin/env python3
"""
File:           point.py
Description:    Class definition for Point used by NetworkX to create figure.
                Inspiration taken from a colleague's plotting code:
                    https://github.com/yingying-m/protein-folding.
"""


class Point:
    """Point definition used for the front pages figure with NetworkX."""

    def __init__(self, x, y):
        self.x = x
        self.y = y

    # Print the coordinates of some point
    def __repr__(self):
        return "".join(["(", str(self.x), ",", str(self.y), ")"])

    # Find the (un)occupied positions adjacent to some point
    def find_positions(self, grid, relative_position_list, is_occupied):
        position_list = []
        for relative_position in relative_position_list:
            adjacent_position = Point(
                self.x + relative_position[0], self.y + relative_position[1]
            )
            if (
                adjacent_position.x > (len(grid) - 1)
                or adjacent_position.x < 0
                or adjacent_position.y > (len(grid[len(grid) - 1]) - 1)
                or adjacent_position.y < 0
            ):
                continue
            if is_occupied:
                if grid[adjacent_position.y][adjacent_position.x] == " ":
                    continue
            else:
                if grid[adjacent_position.y][adjacent_position.x] != " ":
                    continue
            position_list.append(adjacent_position)
        return position_list
