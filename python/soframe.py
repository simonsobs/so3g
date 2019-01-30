from spt3g import core

class SOFrame(core.G3Frame):
	"""Convenience wrapper for G3Frame. It's purpose is to allow storage for non-C++-based
	objects in a G3Frame, such as numpy arrays or enmaps. Classes must still be defined on
	the C++ side, but they don't have to be the same as the corresponding classes on the python
	side. This is handled by defining a set of converters that transparently convert between
	the python and C++ representation. Types not handled by converters are handled as before.
	
	Conversion is controlled via the static getitem_converters and setitem_converters
	members of the SOFrame class. Here's an example of how to register a new type:

	SOFrame.setitem_converters[np.ndarray] = lambda a: G3Ndarray(a)
	SOFrame.getitem_converters[G3Ndarray]  = lambda a: a.get_array()

	With this set, one can now assign the given type (in this case a numpy array)
	to a frame like this:

	frame = SOFrame()
	frame["data"] = np.arange(10.)
	print(frame["data"])
	"""
	getitem_converters = {}
	setitem_converters = {}
	def __getitem__(self, key):
		val = core.G3Frame.__getitem__(self, key)
		try: val = SOFrame.getitem_converters[type(val)](val)
		except KeyError: pass
		return val
	def __setitem__(self, key, val):
		try: val = SOFrame.setitem_converters[type(val)](val)
		except KeyError: pass
		core.G3Frame.__setitem__(self, key, val)
