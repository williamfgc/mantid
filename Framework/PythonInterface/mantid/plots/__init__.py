# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantid package
#
#
"""
Functionality for unpacking mantid objects for plotting with matplotlib.
"""

# This file should be left free of PyQt imports to allow quick importing
# of the main package.
from __future__ import (absolute_import, division, print_function)

from collections import Iterable

import mantid.kernel
import mantid.plots.plotfunctions
import mantid.plots.plotfunctions3D
from matplotlib.axes import Axes
from matplotlib.container import Container
from matplotlib.projections import register_projection

try:
    from mpl_toolkits.mplot3d.axes3d import Axes3D
except ImportError:
    # Special case to handle issues with importing mpl_toolkits
    #
    # Matplotlib adds a *nspkg.pth file to the user site packages directory.
    # When that file is processed a fake built-in mpl_toolkits is imported
    # which forces the site packages version to take precidence over our
    # local copy regardless of python sys path settings.
    #
    # Work around by removing the fake built-in module from sys modules,
    # then forcing python to search the path as expected.
    #
    # This is mostly but not necessarily limited to being an issue on OSX
    # where there are multiple versions of matplotlib installed across the
    # system.
    import sys
    del sys.modules['mpl_toolkits']
    from mpl_toolkits.mplot3d.axes3d import Axes3D


def plot_decorator(func):
    def wrapper(self, *args, **kwargs):
        func_value = func(self, *args, **kwargs)
        # Saves saving it on array objects
        if mantid.plots.helperfunctions.validate_args(*args, **kwargs):
            # Fill out kwargs with the values of args
            for index, arg in enumerate(args):
                if index is 0:
                    kwargs["workspaces"] = args[0].name()
                if index is 1:
                    kwargs["spectrum_nums"] = args[1]
                if index is 2:
                    kwargs["wksp_indices"] = args[2]
                if index is 3:
                    kwargs["errors"] = args[3]
                if index is 4:
                    kwargs["overplot"] = arg[4]
                # ignore 5 as no need to save the fig object
                if index is 6:
                    kwargs["plot_kwargs"] = arg[6]
            if hasattr(self, "creation_args"):
                self.creation_args.append(kwargs)
            else:
                self.creation_args = [kwargs]
        return func_value
    return wrapper


class MantidAxes(Axes):
    """
    This class defines the **mantid** projection for 2d plotting. One chooses
    this projection using::

        import matplotlib.pyplot as plt
        from mantid import plots
        fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})

    or::

        import matplotlib.pyplot as plt
        from mantid import plots
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='mantid')

    The mantid projection allows replacing the array objects with mantid workspaces.
    """
    # Required by Axes base class
    name = 'mantid'

    # Enumerators for plotting directions
    HORIZONTAL = BIN = 0
    VERTICAL = SPECTRUM = 1

    # Store information for any workspaces attached to this axes instance
    tracked_workspaces = None

    def __init__(self, *args, **kwargs):
        super(MantidAxes, self).__init__(*args, **kwargs)
        self.tracked_workspaces = dict()

    def track_workspace_artist(self, name, artists):
        """
        Add the given workspace name to the list of workspaces
        displayed on this Axes instance
        :param name: The name of the workspace
        :param artists: A single artist or list of artists held in the container
        :returns: The artists variable as it was passed in.
        """
        if name:
            artist_info = self.tracked_workspaces.setdefault(name, [])
            if isinstance(artists, Iterable) and not isinstance(artists, Container):
                for artist in artists:
                    artist_info.append(artist)
            else:
                artist_info.append(artists)

        return artists

    def remove_workspace_artists(self, name):
        """
        Remove the artists reference by this workspace (if any) and return True
        if the axes is then empty
        :param name: The name of the workspace
        :return: True if the axes is empty, false if artists remain or this workspace is not associated here
        """
        try:
            artist_info = self.tracked_workspaces.pop(name)
        except KeyError:
            return False

        # delete the artists from the figure
        for artist in artist_info:
            artist.remove()
            # Remove doesn't catch removing the container for errorbars etc
            if isinstance(artist, Container):
                try:
                    self.containers.remove(artist)
                except ValueError:
                    pass

        axes_empty = self.is_empty()
        if (not axes_empty) and self.legend_ is not None:
            self.legend()

        return axes_empty

    def is_empty(self):
        """
        Checks the known artist containers to see if anything exists within them
        :return: True if no artists exist, false otherwise
        """
        def _empty(container):
            return len(container) == 0
        return _empty(self.lines) and _empty(self.images) and _empty(self.collections) and _empty(self.containers)

    @plot_decorator
    def plot(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.plot` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.plot(workspace,'rs',specNum=1) #for workspaces
            ax.plot(x,y,'bo')                 #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.plot`.
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.plot(self, *args, **kwargs)
        else:
            return Axes.plot(self, *args, **kwargs)

    def scatter(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.scatter` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.scatter(workspace,'rs',specNum=1) #for workspaces
            ax.scatter(x,y,'bo')                 #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.scatter`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.scatter(self, *args, **kwargs)
        else:
            return Axes.scatter(self, *args, **kwargs)

    def errorbar(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.errorbar` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.errorbar(workspace,'rs',specNum=1) #for workspaces
            ax.errorbar(x,y,yerr,'bo')            #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.errorbar`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.errorbar(self, *args, **kwargs)
        else:
            return Axes.errorbar(self, *args, **kwargs)

    def pcolor(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.pcolor` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.pcolor(workspace) #for workspaces
            ax.pcolor(x,y,C)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.pcolor`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.pcolor(self, *args, **kwargs)
        else:
            return Axes.pcolor(self, *args, **kwargs)

    def pcolorfast(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.pcolorfast` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.pcolorfast(workspace) #for workspaces
            ax.pcolorfast(x,y,C)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.pcolorfast`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.pcolorfast(self, *args, **kwargs)
        else:
            return Axes.pcolorfast(self, *args, **kwargs)

    def pcolormesh(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.pcolormesh` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.pcolormesh(workspace) #for workspaces
            ax.pcolormesh(x,y,C)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.pcolormesh`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.pcolormesh(self, *args, **kwargs)
        else:
            return Axes.pcolormesh(self, *args, **kwargs)

    def imshow(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.imshow` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.imshow(workspace) #for workspaces
            ax.imshow(C)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.imshow`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.imshow(self, *args, **kwargs)
        else:
            return Axes.imshow(self, *args, **kwargs)

    def contour(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.contour` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.contour(workspace) #for workspaces
            ax.contour(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.contour`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.contour(self, *args, **kwargs)
        else:
            return Axes.contour(self, *args, **kwargs)

    @plot_decorator
    def contourf(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.contourf` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.contourf(workspace) #for workspaces
            ax.contourf(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.contourf`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.contourf(self, *args, **kwargs)
        else:
            return Axes.contourf(self, *args, **kwargs)

    def tripcolor(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.tripcolor` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.tripcolor(workspace) #for workspaces
            ax.tripcolor(x,y,C)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.tripcolor`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.tripcolor(self, *args, **kwargs)
        else:
            return Axes.tripcolor(self, *args, **kwargs)

    def tricontour(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.tricontour` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.tricontour(workspace) #for workspaces
            ax.tricontour(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.tricontour`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.tricontour(self, *args, **kwargs)
        else:
            return Axes.tricontour(self, *args, **kwargs)

    def tricontourf(self, *args, **kwargs):
        """
        If the **mantid** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes.tricontourf` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid'})
            ax.tricontourf(workspace) #for workspaces
            ax.tricontourf(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions.tricontourf`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions')
            return mantid.plots.plotfunctions.tricontourf(self, *args, **kwargs)
        else:
            return Axes.tricontourf(self, *args, **kwargs)


class MantidAxes3D(Axes3D):
    """
    This class defines the **mantid3d** projection for 3d plotting. One chooses
    this projection using::

        import matplotlib.pyplot as plt
        from mantid import plots
        fig, ax = plt.subplots(subplot_kw={'projection':'mantid3d'})

    or::

        import matplotlib.pyplot as plt
        from mantid import plots
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='mantid3d')

    The mantid3d projection allows replacing the array objects with mantid workspaces.
    """

    name = 'mantid3d'

    def plot(self, *args, **kwargs):
        """
        If the **mantid3d** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes3D.plot` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid3d'})
            ax.plot(workspace) #for workspaces
            ax.plot(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions3D.plot3D`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions3D')
            return mantid.plots.plotfunctions3D.plot(self, *args, **kwargs)
        else:
            return Axes3D.plot(self, *args, **kwargs)

    def scatter(self, *args, **kwargs):
        """
        If the **mantid3d** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes3D.scatter` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid3d'})
            ax.scatter(workspace) #for workspaces
            ax.scatter(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions3D.scatter`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions3D')
            return mantid.plots.plotfunctions3D.scatter(self, *args, **kwargs)
        else:
            return Axes3D.scatter(self, *args, **kwargs)

    def plot_wireframe(self, *args, **kwargs):
        """
        If the **mantid3d** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes3D.plot_wireframe` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid3d'})
            ax.plot_wireframe(workspace) #for workspaces
            ax.plot_wireframe(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions3D.wireframe`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions3D')
            return mantid.plots.plotfunctions3D.plot_wireframe(self, *args, **kwargs)
        else:
            return Axes3D.plot_wireframe(self, *args, **kwargs)

    def plot_surface(self, *args, **kwargs):
        """
        If the **mantid3d** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes3D.plot_surface` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid3d'})
            ax.plot_surface(workspace) #for workspaces
            ax.plot_surface(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions3D.plot_surface`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions3D')
            return mantid.plots.plotfunctions3D.plot_surface(self, *args, **kwargs)
        else:
            return Axes3D.plot_surface(self, *args, **kwargs)

    def contour(self, *args, **kwargs):
        """
        If the **mantid3d** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes3D.contour` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid3d'})
            ax.contour(workspace) #for workspaces
            ax.contour(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions3D.contour`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions3D')
            return mantid.plots.plotfunctions3D.contour(self, *args, **kwargs)
        else:
            return Axes3D.contour(self, *args, **kwargs)

    def contourf(self, *args, **kwargs):
        """
        If the **mantid3d** projection is chosen, it can be
        used the same as :py:meth:`matplotlib.axes.Axes3D.contourf` for arrays,
        or it can be used to plot :class:`mantid.api.MatrixWorkspace`
        or :class:`mantid.api.IMDHistoWorkspace`. You can have something like::

            import matplotlib.pyplot as plt
            from mantid import plots

            ...

            fig, ax = plt.subplots(subplot_kw={'projection':'mantid3d'})
            ax.contourf(workspace) #for workspaces
            ax.contourf(x,y,z)     #for arrays
            fig.show()

        For keywords related to workspaces, see :func:`mantid.plots.plotfunctions3D.contourf`
        """
        if mantid.plots.helperfunctions.validate_args(*args):
            mantid.kernel.logger.debug('using mantid.plots.plotfunctions3D')
            return mantid.plots.plotfunctions3D.contourf(self, *args, **kwargs)
        else:
            return Axes3D.contourf(self, *args, **kwargs)


register_projection(MantidAxes)
register_projection(MantidAxes3D)
