from .base import Parser
from ..utils.constants import C_LIGHT
import matplotlib.pyplot as plt
import numpy as np

class Bandstructure(Parser):
  """Holds the result of a bandstructure simulation"""
  def show(self, title="Bandstructure"):
    """Plot the bandstructure as a pcolormesh

    Parameters
    ----------
    title : str
      The title of the plot
    """
    krange = self.meta[0]
    freqrange = self.meta[1]
    a = self.meta[2]

    res2 = self.data.copy()
    for i in range(res2.shape[0]):
      res2[i,:] = res2[i,:] / np.max(res2[i,:])
    res2 = np.transpose(res2)

    fig, ax = plt.subplots()
    c = ax.pcolormesh(krange, freqrange, res2)
    plt.colorbar(c, ax=[ax], location="left")

    lightlinex = np.array([krange[0], krange[-1]])
    if krange[0] * C_LIGHT / a < freqrange[0]:
      lightlinex[0] = a*freqrange[0]/C_LIGHT
    if krange[-1] * C_LIGHT / a > freqrange[-1]:
      lightlinex[1] = a*freqrange[-1]/C_LIGHT
    ax.plot(lightlinex, lightlinex * C_LIGHT / a, color="black")

    ax.set_xlabel("K")
    ax.set_ylabel("Frequency")
    ax.get_yaxis().set_ticks_position("right")
    ax.get_yaxis().set_label_position("right")
    plt.title(title)

    plt.show()
  
  def __repr__(self):
    return "Bandstructure((%.6e, %.6e), (%.6e, %.6e), %.6e)" % (self.meta[0][0], self.meta[0][-1], self.meta[1][0], self.meta[1][-1], self.meta[2])

  def save(self, fpath):
    # TODO: Implement
    pass

def correctDiscontinuity(field):
    '''When using symetries, for some reason, if the max value is a unic value in the matrix,
        is a unusual high value. This function removes this value if it is the case'''
      
    bm = field == field.max()
    if bm.sum() == 1:
      print('Removing possible discontinuity due symmetries')
      posMax = np.unravel_index(np.argmax(field),field.shape)
      field[posMax[0],posMax[1]] = field[posMax[0]+1,posMax[1]+1]

class EField(Parser):
  """Holds an electric field profile result from a cavity resonance simulation"""
  def show(self, title=None, ncontours=1,show = True,cmapFields = 'seismic'):
    """Plots the profile in 4 graphs, one for each of the x, y, and z components of the
    electric field, and one for the magnitude of the field. This also overlays a contour
    graph of the index of refraction of the cavity to show where in the geometry the field is 
    concentrated

    Parameters
    ----------
    title : str or None
      If provided, override the default title for the plot
    ncontours : int
      Number of contour lines used in the index of refraction contour plot
    """
    xlabel = self.meta[0]
    ylabel = self.meta[1]

    if title is None:
      title = "Electric field " + (xlabel + ylabel).upper() + " density profile"
    
    Ex, Ey, Ez, x, y, index = self.data
    E = Ex+ Ey + Ez
    fig, axs = plt.subplots(2, 2)
    for cname, data, ax in [
        ("Ex", Ex, axs[0,0]), ("Ey", Ey, axs[0,1]), ("Ez", Ez, axs[1,0]), ("E", E, axs[1,1])
      ]:
      correctDiscontinuity(data)
      c = ax.pcolormesh(x, y, np.transpose(data),cmap = cmapFields)
      ax.contourf(x, y, np.transpose(index), ncontours, cmap="binary",alpha = 0.2)
##      ax.contourf(x, y, np.transpose(index), ncontours, colors="black")
      plt.colorbar(c, ax=[ax], location="left")

      ax.set_xlabel(xlabel)
      ax.set_ylabel(ylabel)
      ax.set_aspect("equal")
      ax.get_yaxis().set_ticks_position("right")
      ax.get_yaxis().set_label_position("right")
      ax.set_title(title + " (" + cname + ")")
    if show == True:
      plt.show()
    else:
      return fig
    
  def __repr__(self):
    return "EField(%s, %s)" % (self.meta[0], self.meta[1])

  def max_loc(self):
    Ex, Ey, Ez, x, y, index = self.data
    E = Ex + Ey + Ez
    print(np.shape(Ex))
    print(np.shape(Ey))
    print(np.shape(Ez))
    print(np.shape(E))
    print(np.shape(x))
    print(np.shape(y))
    max_index = np.unravel_index(np.argmax(E),E.shape)
    print(max_index)
    return (x[max_index[0]],y[max_index[1]])
    

  def save(self, fpath):
    self.show(show = False).savefig(fpath)
    pass
    
class Quasipotential(Parser):
  """Holds the result of a quasipotential simulation"""
  def show(self, title="Quasipotential"):
    """Plot the quasipotential

    Parameters
    ----------
    title : str
      Title of the plot
    """
    plt.plot([0, len(self.data) - 1], [0, 0], color="black")
    plt.title(title)
    plt.ylabel("Potential (THz)")
    plt.xlabel("Cell index")
    plt.xticks(list(range(len(self.data))))
    plt.plot(self.data)
    plt.show()
  
  def __repr__(self):
    return "Quasipotential(%d)" % len(self.data)

  def save(self, fpath):
    # TODO: Implement
    pass
