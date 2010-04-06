"""
Simulated annealing with restart applied to the traveling salesman problem.

The data needs to be in cvs format and the coordinates of the cities are
longitude and latitude. Thus, any gps data taken on the internet can
serve as input.

Once the annealing is done, the cities are plotted according to their
coordinates using the matplotlib package. The cities are then linked as to
show the order of visit in the final solution. Also, the evolution of the
tested and shortest distances is plotted.
"""
__docformat__ = "restructuredtext en" 

## Copyright (c) 2010 Emmanuel Goossaert 
##
## This file is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.


import csv
import sys
import math
import random
import time

import matplotlib.pyplot as plt


def read_cities(filename):
    """
    Read city data from a csv file.
    """
    reader = csv.reader(open(filename, "rb")) # may raise IOError
    rows = [line for line in reader]
    cities = [City(r[2], r[3], float(r[0]), float(r[1])) for r in rows[1:]]
    return cities



class City:
    """
    Store information regarding a city, including name and gps coordinates.
    """

    def __init__(self, name='', description='', latitude=0, longitude=0):
        self.name = name
        self.description = description
        self.latitude = latitude
        self.longitude = longitude

    def __str__(self):
        return self.name + ' ' + self.description + ' ' + str(self.latitude) + ' ' + str(self.longitude)

    def __repr__(self):
        return self.__str__()

    def distance_to_city_in_km(self, city):
        """Distance to another city using Haversine formula"""
        lat = math.radians(self.latitude - city.latitude)
        long = math.radians(self.longitude - city.longitude)
        a = math.pow(math.sin(lat/2), 2) \
           + math.cos(math.radians(self.latitude)) * math.cos(math.radians(city.latitude)) * pow(math.sin(long/2), 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

        radius_earth = 6378.7 # in kilometers
        return radius_earth * c


def total_distance_in_km(cities):
    distances = [cities[index].distance_to_city_in_km(cities[(index + 1) % len(cities)]) for index in range(len(cities))]
    return sum(distances)


def plot_cities(cities, figure_id):
    """Plot the cities on a plan."""
    fig_map = plt.figure(figure_id)
    ax_map = fig_map.add_subplot(111)

    cities_x = [city.longitude for city in cities + [cities[0]]]
    cities_y = [city.latitude for city in cities + [cities[0]]]

    link = '-'
    ax_map.plot(cities_x, cities_y, 'go' + link)
    ax_map.grid()

    spacing = math.fabs(min(cities_x) - max(cities_x)) * .1
    ax_map.set_xlim(min(cities_x) - spacing, max(cities_x) + spacing * 3)
    ax_map.set_ylim(min(cities_y) - spacing, max(cities_y) + spacing)

    for index, city in enumerate(cities):
        ax_map.text(city.longitude,
                    city.latitude,
                    '%d: %s' % (index + 1, city.name),
                    withdash = True,
                    dashdirection = 1,
                    dashlength = 30,
                    rotation = 0,
                    dashrotation = 15,
                    dashpush = 10)
    return ax_map


def plot_distances(distances_current, figure_id, distances_best, ids_restart, nb_cities, nb_iterations):
    """Plot the evolution of the distance metrics."""
    # plot distances
    fig_distances = plt.figure(figure_id)
    ax_distances = fig_distances.add_subplot(111)
    line_current = ax_distances.plot(distances_current, linewidth=1)
    line_best = ax_distances.plot(distances_best, 'r', linewidth=2)
    ax_distances.set_title('Distance evolution for %d cities on %d iteration(s)' % (nb_cities, nb_iterations))

    # plot restart steps
    y_min = min(distances_current)
    y_max = max(distances_current)
    line_restart = None

    for step in ids_restart[:-1]:
        line_restart = ax_distances.plot([step, step], [y_min, y_max], 'g', linewidth=2)

    ax_distances.set_xlabel('Steps (with systematic sampling at 1%)')
    ax_distances.set_ylabel('Distance (km)')

    index_legend = 3 if len(ids_restart) > 1 else 2
    plt.legend( (line_current, line_best, line_restart)[:index_legend],
                ('Tested distance', 'Shortest distance', 'Restarting annealing process')[:index_legend],
                loc='upper right' )
    

def annealing(cities, temperature_begin=100000000.0, temperature_end=.1, cooling_factor=.95, nb_iterations=1):
    """
    Simulated annealing function, implemented with acceptance probability from
    by Kirkpatrick et al., and with restart.

    distance_best:    best solution encountered so far
    distance_current: solution used in the current simulation
    distance_new:     solution computed from the random changes to current
    """

    cities_best = cities[:]
    distance_best = total_distance_in_km(cities_best)

    distances_current = []
    distances_best = []
    ids_restart = []

    try:
        for iteration in range(nb_iterations):
            # the search is restarted at every iteration from
            # the best know solution
            temperature = temperature_begin
            cities_current = cities_best[:]
            distance_current = distance_best
            cities_new = cities_best[:]

            step = 0
            while temperature > temperature_end:
                # swap two random cities, but never touch the first city
                # -- could be optimized by recomputing only the changed distances
                index = random.sample(range(len(cities_new) - 1), 2)
                index[0] += 1
                index[1] += 1
                cities_new[index[0]], cities_new[index[1]] = cities_new[index[1]], cities_new[index[0]]

                # compute the new distance
                distance_new = total_distance_in_km(cities_new)

                # acceptance probability by Kirkpatrick et al.
                diff = distance_new - distance_current
                if diff < 0 or  math.exp( -diff / temperature ) > random.random():
                    cities_current = cities_new[:]
                    distance_current = distance_new

                # update the best if current solution is better
                # not part of the annealing itself, just used for the restart
                if distance_current < distance_best:
                    cities_best = cities_current[:]
                    distance_best = distance_current
                    step = 0 # reset step to make sure the point is saved

                if step % 100 == 0:
                    # systematic sampling: one point every 100th
                    distances_current.append(distance_current)
                    distances_best.append(distance_best)
                temperature = temperature * cooling_factor
                step = step + 1

            ids_restart.append(len(distances_current))

    except KeyboardInterrupt, e:
        print "Interrupted on user demand."
        print 'performed iterations: %d' % iteration

    return cities_best, distances_current, distances_best, ids_restart


def display_usage():
    print 'usage: %s input [nb_ite] [nb_cities] [c_factor] [t_start] [t_end]' % sys.argv[0]
    print '  input: input CSV file containing the city coordinates'
    print '  nb_ite: number of iterations in the restart process'
    print '  nb_cities: number of cities read from the input file (n first lines)'
    print '  c_factor: cooling factor, float number in (0,1)'
    print '  t_start: initial temperature'
    print '  t_end: temperature at which the process will be stopped'
    print '         must be smaller than t_start'


if __name__ == '__main__':

    if len(sys.argv) < 2:
        display_usage()
        sys.exit(0)

    # Initialize random number generator
    # during development, keep constant seed value: 42
    random.seed(42) 

    input = sys.argv[1]
    nb_iterations     = int(sys.argv[2])   if len(sys.argv) > 2 else 1
    nb_cities         = int(sys.argv[3])   if len(sys.argv) > 3 else -1
    cooling_factor    = float(sys.argv[4]) if len(sys.argv) > 4 else .99
    temperature_start = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0e+300
    temperature_end   = float(sys.argv[6]) if len(sys.argv) > 6 else .1

    time_begin = time.time()
    cities = read_cities(input)

    nb_cities = len(cities) if nb_cities <= 0 else nb_cities

    cities = cities[:nb_cities]
    (cities_new, distances_current, distances_best, ids_restart) = annealing(cities, temperature_start, temperature_end, cooling_factor, nb_iterations)
    time_end = time.time()

    distance_begin = total_distance_in_km(cities)
    distance_end = total_distance_in_km(cities_new)
    print 'Improvement:          %8.0f %%'  % (100 * (distance_begin - distance_end) / distance_begin)
    print 'Time:                 %8.0f sec' % (time_end - time_begin)
    print 'Initial distance:     %8.0f km'  % distance_begin
    print 'Final distance:       %8.0f km'  % distance_end

    ax_map = plot_cities(cities, 1)
    ax_map.set_title('Initial distance: %.0f km for %d cities' % (distance_begin, len(cities)))

    if nb_iterations:
        ax_map = plot_cities(cities_new, 2)
        ax_map.set_title('Final distance: %.0f km for %d cities on %d iteration(s)' % (distance_end, len(cities), nb_iterations))
        plot_distances(distances_current, 3, distances_best, ids_restart, len(cities), nb_iterations)

    plt.show()
