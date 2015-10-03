#!/usr/bin/python
import logging
import csv
import numpy as np

from copy import deepcopy


class BaseBounds(object):
    """A base class for declaring bounds on things."""

    def Copy(self):
        """Returns a (deep) copy of self."""
        raise NotImplementedError

    def GetRange(self):
        """Returns a 2-tuple of the concentration range."""
        return None

    def GetLowerBound(self, key):
        """Get the lower bound for this key.
        
        Args:
            key: a string key.
        """
        raise NotImplementedError
        
    def GetUpperBound(self, key):
        """Get the upper bound for this key.
        
        Args:
            key: a string key.
        """
        raise NotImplementedError

    def GetBounds(self, keys):
        """Get the bounds for a set of keys in order.
        
        Args:
            keys: an iterable of keys.
        
        Returns:
            A two-tuple (lower_bounds, upper_bounds) where both
            items are Numpy arrays of dimensions 1xlen(keys)
        """
        lower_bounds = np.array([self.GetLowerBound(key) for key in keys])
        upper_bounds = np.array([self.GetUpperBound(key) for key in keys])
        lower_bounds = np.matrix(lower_bounds.reshape(1, len(lower_bounds)))
        upper_bounds = np.matrix(upper_bounds.reshape(1, len(upper_bounds))) 
        return lower_bounds, upper_bounds
    
    def GetBoundsWithDefault(self, keys, default):
        """Returns the default value for each key unless it is invalid.
        
        Returns:
            A row vector (Numpy array) with shape 1xlen(keys).
        """
        res = []
        for key in keys:
            ub = self.GetUpperBound(key)
            lb = self.GetLowerBound(key)
            
            if default < lb:
                res.append(lb)
                continue
            
            if default > ub:
                res.append(ub)
                continue
            
            res.append(default)
        
        a = np.array(res)
        return a.reshape(1, len(res))
        
    def GetLnBounds(self, keys):
        """Get the bounds for a set of keys in order.
        
        Args:
            keys: an iterable of keys.
        
        Returns:
            A two-tuple (lower_bounds, upper_bounds) where both
            items are Numpy arrays of dimensions 1xlen(keys)
        """
        lb, ub = self.GetBounds(keys)
        return np.log(lb), np.log(ub)
    
    def SetBounds(self, key, lb, ub):
        """Set bounds for a specific key.
        
        Args:
            key: the key for the bounds.
            lb: the lower bound value.
            ub: the upper bound value.
        """
        assert lb <= ub
        self.lower_bounds[key] = lb
        self.upper_bounds[key] = ub
    
    def AddOldStyleBounds(self, bounds_dict):
        """Add bounds from the old-style bounds dictionary.
        
        Args:
            bounds_dict: a dictionary mapping keys to (lb, ub) tuples.
        """
        for key, bounds in bounds_dict.iteritems():
            lb, ub = bounds
            self.lower_bounds[key] = lb
            self.upper_bounds[key] = ub
    
    def GetOldStyleBounds(self, keys):
        lower_bounds, upper_bounds = self.GetBounds(keys)
        cid2bounds = {}
        for i, cid in enumerate(keys):
            cid2bounds[cid] = (lower_bounds[0, i], upper_bounds[0, i])
        return cid2bounds


class ExplicitBounds(BaseBounds):
    """Contains upper and lower bounds for various keys."""

    def __init__(self, lower_bounds=None, upper_bounds=None):
        """Initialize the ExplicitBounds object.
        
        Must provide upper and lower bounds for all compounds.
        
        Args:
            lower_bounds: a dictionary mapping strings to float lower bounds.
            upper_bounds: a dictionary mapping strings to float upper bounds.
        """        
        self.lower_bounds = lower_bounds or {}
        self.upper_bounds = upper_bounds or {}
        
        # Must have the same keys for both
        lb_keys = set(self.lower_bounds.keys())
        ub_keys = set(self.upper_bounds.keys())
        assert lb_keys == ub_keys

        self.c_range = (np.min(self.lower_bounds), np.max(self.upper_bounds))

    def Copy(self):
        """Returns a deep copy of self."""
        new_lb = deepcopy(self.lower_bounds)
        new_ub = deepcopy(self.upper_bounds)
        return ExplicitBounds(new_lb, new_ub)

    def GetRange(self):
        """Returns a 2-tuple of the concentration range."""
        return self.c_range

    def GetLowerBound(self, key):
        """Get the lower bound for this key.
        
        Args:
            key: a string key.
        """
        if key not in self.lower_bounds:
            raise KeyError('Unknown key %s' % key)
        
        return self.lower_bounds[key]
    
    def GetUpperBound(self, key):
        """Get the upper bound for this key.
        
        Args:
            key: a string key.
        """
        if key not in self.upper_bounds:
            raise KeyError('Unknown key %s' % key)
        
        return self.upper_bounds[key]
    
    

class Bounds(BaseBounds):
    """Contains upper and lower bounds for various keys."""
    
    def __init__(self,
                 lower_bounds=None,
                 upper_bounds=None,
                 default_lb=None,
                 default_ub=None):
        """Initialize the bounds object.
        
        Args:
            lower_bounds: a dictionary mapping strings to float lower bounds.
            upper_bounds: a dictionary mapping strings to float upper bounds.
            default_lb: the default lower bound to return.
            default_lb: the default upper bound to return.
        """
        self.lower_bounds = lower_bounds or {}
        self.upper_bounds = upper_bounds or {}
        self.default_lb = default_lb
        self.default_ub = default_ub
        
        self.c_range = (self.default_lb, self.default_ub)

    @classmethod
    def from_csv_file(cls, f, default_lb=None, default_ub=None):
        lbs = {}
        ubs = {}
        for row in csv.DictReader(f):
            cid = row['cid']
            lb, ub = row.get('c_min'), row.get('c_max')
            if lb.strip():
                lb = float(lb)
            else:
                lb = None
            if ub.strip():
                ub = float(ub)
            else:
                ub = None

            lbs[cid] = lb
            ubs[cid] = ub
        return Bounds(lbs, ubs, default_lb, default_ub)

    @classmethod
    def from_csv_filename(cls, fname, default_lb=None, default_ub=None):
        with open(fname) as f:
            return cls.from_csv_file(
                f, default_lb=default_lb, default_ub=default_ub)

    def Copy(self):
        """Returns a deep copy of self."""
        new_lb = deepcopy(self.lower_bounds)
        new_ub = deepcopy(self.upper_bounds)
        return Bounds(new_lb, new_ub,
                      self.default_lb,
                      self.default_ub)
        
    def GetRange(self):
        """Returns a 2-tuple of the concentration range."""
        return self.c_range

    def GetLowerBound(self, key):
        """Get the lower bound for this key.
        
        Args:
            key: a string key.
        """
        val = self.lower_bounds.get(key) or self.default_lb
        return val
    
    def GetUpperBound(self, key):
        """Get the upper bound for this key.
        
        Args:
            key: a string key.
        """
        val = self.upper_bounds.get(key) or self.default_ub
        return val
    