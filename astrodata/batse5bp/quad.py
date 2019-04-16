"""
Generic basic and composite quadrature rule objects.

Created 2012-10-28 by Tom Loredo; stolen from the Inference package
2019:  Converted to Python 3
"""

from numpy import array, concatenate, searchsorted


class Quad:
    """
    A simple quadrature rule
    """

    def __init__(self, *args):
        """
        Define a quadrature rule from its nodes and weights.

        For closed rules (where the nodes include the boundaries), the
        signature may be:

            QuadRule(nodes, wts)

        where `nodes` and `wts` are arrays containing the (ordered)
        nodes and weights for the rule.

        For open rules (and optionally for closed rules), the signature is

            QuadRule(l, u, nodes, wts)

        where `a` and `b` are the integration limits, and `nodes` and `wts`
        are arrays containing the (ordered) nodes and weights for the rule.
        """
        if len(args) == 2:
            self.type = 'closed'
            self.nodes = args[0]
            self.wts = args[1]
            self.l = self.nodes[0]
            self.u = self.nodes[-1]
        elif len(args) == 4:
            self.l = args[0]
            self.u = args[1]
            self.nodes = args[2]
            self.wts = args[3]
            if self.l == self.nodes[0] and self.u == self.nodes[-1]:
                self.type = 'closed'
            else:
                self.type = 'open'
        else:
            raise ValueError('Specify (nodes, wts) or (a, b, nodes, wts)!')
        self.npts = len(self.nodes)
        if len(self.wts) != self.npts:
            raise ValueError('Mismatch in length of nodes and wts!')

    def quad(self, f):
        """
        Evaluate the quadrature given a function or array of function values.

        If f is callable, find the array of values f(self.nodes) and
        evaluate the quadrature.  Otherwise, f must be a vector of function
        values at the nodes, used to evaluate the quadrature.
        """
        if callable(f):
            fvals = f(self.nodes)
        elif len(f) == self.npts:
            fvals = f
        else:
            raise ValueError('Argument must be callable or array of values!')
        return sum(self.wts*fvals)


class CompositeQuad:
    """
    Composite quadrature rule built over contiguous intervals
    """

    @staticmethod
    def isquad(obj):
        """
        Return True if obj has a Quad interface.

        Note that this does not check for a quad_range method; this will
        not be needed by some CompositeQuad instances.
        """
        if isinstance(obj, Quad):
            return True
        try:
            assert hasattr(obj, 'l')
            assert hasattr(obj, 'u')
            assert hasattr(obj, 'npts')
            assert hasattr(obj, 'nodes')
            assert hasattr(obj, 'wts')
            return True
        except AssertionError:
            return False

    def __init__(self, *args, **kwds):
        """
        Define a composite quadrature rule from a collection of rules.

        Each argument should describe a basic quadrature rule.  There are
        three possible formats:

          * The argument may be a QuadRule instance.

          * For a closed rule, the argument may be a 2-tuple of the
            form (nodes, wts), where `nodes` and `wts` are arrays of
            values for the nodes and weights of the rule.

          * For an open or closed rule, the argument may be a 4-tuple of the
            form (a, b, nodes, wts), where `a` and `b` are the limits of
            integration, and `nodes` and `wts` are as above.

        The rules must be contiguous.

        If there is a `factor` keyword argument, it should be a float considered
        to be a constant factor taken out of the integrand; it will be applied
        at the end of quadrature calculations.

        Quadrature over part of the region spanned by the provided rules is
        supported.  This requires the rules to have methods that integrate
        over part of their range, and that return nodes and weights for part
        of their range.  By default, the method name to be used for quadrature
        is 'quad_range', but a 'range_method' keyword argument may be provided,
        specifying the name of the range-altering quadrature method.  The
        method name used for getting new nodes and weights can be similarly
        specified with a 'nw_method' keyword argument, with the default name
        being 'quad_range_nodes_wts'.
        """
        # Collect all the rules as QuadRule instances.
        self.rules = []
        last = None
        for i, arg in enumerate(args):
            if CompositeQuad.isquad(arg):
                if last:
                    if last.u != arg.l:
                        raise ValueError(
                            'Rule %i not contiguous with previous!' % (i+1))
                self.rules.append(arg)
                last = arg
            else:
                rule = Quad(*arg)
                self.rules.append(rule)
                last = rule
        self.n_rules = len(self.rules)

        # We can compute the quadrature by running through the rules and
        # computing each one separately, but to accelerate the computation we
        # make an array of all the distinct nodes and their associated weights,
        # summing weights of repeated nodes.
        prev = self.rules[0]
        nodes = [prev.nodes]  # we'll later concatenate these
        wts = [prev.wts.copy()]    # copy for when boundary wts change
        starts = [0]  # starts[i] = node index for 1st node from rule i
        self.npts = prev.npts
        # To support the ability to alter the integration range, we'll also keep
        # track of the ranges of the rules, rule ids for each node, and the
        # weight contributions for duplicated nodes.
        bounds = [prev.l]
        node_ids = [0]*prev.npts  # gives id of 1st rule containing the node
        node_dup = [False]*prev.npts  # will be tuple of wts if node is a dup
        for i, rule in enumerate(self.rules[1:]):  # note i = rule # - 1
            bounds.append(rule.l)
            if prev.nodes[-1] != rule.nodes[0]:  # no boundary node duplication
                nodes.append(rule.nodes)
                wts.append(rule.wts.copy())
                starts.append(self.npts)
                self.npts += rule.npts
                node_ids.extend([i+1]*rule.npts)
                node_dup.extend([False]*rule.npts)
            else:  # boundary node duplication
                new_pts = rule.npts - 1
                nodes.append(rule.nodes[1:])
                node_dup[-1] = (wts[-1][-1], rule.wts[0])  # contributions to wt
                wts[-1][-1] += rule.wts[0]
                wts.append(rule.wts[1:].copy())  # this has to be a copy, too
                starts.append(self.npts-1)
                self.npts += new_pts
                node_ids.extend([i+1]*new_pts)
                node_dup.extend([False]*new_pts)
            prev = rule
        self.bounds = array(bounds + [prev.u])
        self.starts = array(starts)
        self.node_dup = node_dup
        self.nodes = concatenate(nodes)
        self.wts = concatenate(wts)

        # Lower & upper integration bounds:
        self.l = self.rules[0].l
        self.u = self.rules[-1].u

        if 'factor' in kwds:
            self.factor = kwds['factor']
        else:
            self.factor = 1.
        # TODO:  Is there a use case for having separate factors for each rule?
        # self.factors = None

        if 'range_method' in kwds:
            self.rm_name = kwds['range_method']
        else:
            self.rm_name = 'quad_range'
        if 'nw_method' in kwds:
            self.nwm_name = kwds['nw_method']
        else:
            self.nwm_name = 'quad_range_nodes_wts'

    def quad(self, f, call_rules=False):
        """
        Evaluate the quadrature using the callable f or an array of values
        of the integrand at the nodes.

        If call_rules is True, the quadrature is calculated via calls to the
        component rules.  If False, it is calculated using the collected
        distinct nodes and accumulated weights.  These should give exactly
        the same result; the latter (default) approach avoids some overhead
        that can be significant for simple integrands.
        """
        # TODO:  If f() doesn't broadcast, iterate over nodes.
        if callable(f):
            fvals = f(self.nodes)
        elif len(f) == self.npts:
            fvals = f
        else:
            raise ValueError('Argument must be callable or array of values!')
        self.fvals = fvals
        if call_rules:
            # Go through the rules, passing the function evaluations.
            result = 0.
            for rule, start in zip(self.rules, self.starts):
                result += rule.quad(self.ivals[start:start+rule.npts])
            return self.factor*result
        else:
            return self.factor*sum(self.wts*fvals)

    def range_nodes_wts(self, l=None, u=None):
        """
        Return nodes and weights for quadrature over the range [l,u].  If
        either `l` or `u` is None, do not alter that integration bound.
        """
        if l is not None:
            if l < self.l:
                raise ValueError('Lower limit out of range!')
        if u is not None:
            if u > self.u:
                raise ValueError('Upper limit out of range!')

        # Use a subset of the collected nodes and weights, and then construct
        # edge rules for any leftover parts.

        # First, find the ranges for the various parts:  the range spanned by
        # existing nodes, and ranges at each end spanning part of a rule.

        # Note that to find the subset of existing nodes to use, we search over
        # the rule bounds, not the node locations (they will lie inside the
        # bounds for open rules).

        # r_l will point to the left-closed rule containing l.
        print('full range, bounds:', self.l, self.u)
        print(self.bounds)
        if l is None:
            r_l = 0  # lowermost rule to include
            lower = None
        elif l == self.l:
            r_l = 0
            lower = None
        else:
            # Searchsorted gives an insertion index, for the largest bound <= l
            r_l = searchsorted(self.bounds, l)
            r_l -= 1  # since counting intervals, not elements; 0 case handled above
            lower = (l, self.rules[r_l].u)
            if l == lower[1]:  # if l @ right boundary, use the next rule
                lower = None
                r_l += 1
        # r_u will point to the right-closed rule containing u.
        if u is None:
            r_u = self.n_rules - 1  # uppermost rule index
            upper = None
        else:
            r_u = searchsorted(self.bounds, u)  # index of largest bound <= u
            r_u -= 1  # since counting intervals, not elements
            upper = (self.rules[r_u].l, u)
            if u == self.rules[r_u].u:  # if u @ right boundary, use the full rule
                upper = None
        print('l, r_l:', l, r_l, self.rules[r_l].l, self.rules[r_l].u)
        print('u, r_u:', u, r_u, self.rules[r_u].l, self.rules[r_u].u)
        print('lo, up:', lower, upper)

        # Now collect the nodes and weights.
        if r_u == r_l:  # [l,u] lies within 1 rule -> get new nodes, weights
            nw_method = getattr(self.rules[r_l], self.nwm_name)
            if lower:
                a = lower[0]
            else:
                a = self.rules[r_l].l
            if upper:
                b = upper[1]
            else:
                b = self.rules[r_l].u
            return nw_method(a, b)

        else:  # there are full and partial rules; collect & concatenate
            nlist = []
            wlist = []

            # First get the part covered by full rules (if any).
            r_full_l = r_l
            if lower:  # next rule is the lowest possible full one
                r_full_l += 1
            r_full_u = r_u
            if upper:  # prev rule is the highest possible full one
                r_full_u -= 1
            if r_full_u < r_full_l:  # no full rules; [l,u] straddles a boundary
                pass
            else:  # do the full rules
                n_l = self.starts[r_full_l]
                n_u = self.starts[r_full_u] + self.rules[r_full_u].npts - 1
                nodes = self.nodes[n_l:n_u+1]
                if self.node_dup[n_l] or self.node_dup[n_u]:
                    wts = self.wts[n_l:n_u+1].copy()
                    if self.node_dup[n_l]:
                        wts[0] = self.node_dup[n_l][1]
                    if self.node_dup[n_u]:
                        wts[-1] = self.node_dup[n_u][0]
                else:  # don't copy if we don't have to
                    wts = self.wts[n_l:n_u+1]
                # TODO:  If f() doesn't broadcast, iterate over nodes.
                nlist.append(nodes)
                wlist.append(wts)

            # Add partial intervals at each end.
            if lower:
                nw_method = getattr(self.rules[r_l], self.nwm_name)
                nodes, wts = nw_method(lower[0], lower[1])
                print('Lower partial rule nodes, wts:')
                print(nodes)
                print(wts)
                nlist.insert(0, nodes)
                wlist.insert(0, wts)
            if upper:
                nw_method = getattr(self.rules[r_u], self.nwm_name)
                nodes, wts = nw_method(upper[0], upper[1])
                nlist.append(nodes)
                wlist.append(wts)

            return concatenate(nlist), concatenate(wlist)

    def quad_range(self, f, l=None, u=None):
        nodes, wts = self.range_nodes_wts(l, u)
        return self.factor*sum(f(nodes)*wts)

    # *** This is an older version, prior to refactoring out collection of
    # nodes and weights.  Can probably remove; the new version appears to be
    # slightly faster.
    def quad_range0(self, f, l=None, u=None):
        """
        Evaluate the quadrature over [l,u] using the callable f().  If
        either `l` or `u` is None, do not alter that integration bound.
        """
        if l is not None:
            if l < self.l:
                raise ValueError('Lower limit out of range!')
        if u is not None:
            if u > self.u:
                raise ValueError('Upper limit out of range!')
        # Use a subset of the collected nodes and weights, and then apply
        # edge rules for any leftover parts.
        # Note that to find the subset of nodes to use, we search over the
        # rule bounds, not the node locations (they will lie inside the bounds
        # for open rules).
        # r_l will point to the left-closed rule containing l.
        if l is None:
            r_l = 0  # lowermost rule to include
            lower = None
        elif l == self.l:
            r_l = 0
            lower = None
        else:
            # Searchsorted gives an insertion index, for the largest bound <= l
            r_l = searchsorted(self.bounds, l)
            r_l -= 1  # since counting intervals, not elements; 0 case handled above
            lower = (l, self.rules[r_l].u)
            if l == lower[1]:  # if l @ right boundary, use the next rule
                lower = None
                r_l += 1
        # r_u will point to the right-closed rule containing u.
        if u is None:
            r_u = self.n_rules - 1  # uppermost rule index
            upper = None
        else:
            r_u = searchsorted(self.bounds, u)  # index of largest bound <= u
            r_u -= 1  # since counting intervals, not elements
            upper = (self.rules[r_u].l, u)
            if u == self.rules[r_u].u:  # if u @ right boundary, use the full rule
                upper = None

        if r_u == r_l:  # [l,u] lies within 1 rule
            range_method = getattr(self.rules[r_l], self.rm_name)
            if lower:
                a = lower[0]
            else:
                a = self.rules[r_l].l
            if upper:
                b = upper[1]
            else:
                b = self.rules[r_l].u
            result = range_method(f, a, b)
        else:  # there are full and partial rules
            # First do the part covered by full rules (if any).
            r_full_l = r_l
            if lower:  # next rule is the lowest possible full one
                r_full_l += 1
            r_full_u = r_u
            if upper:  # prev rule is the highest possible full one
                r_full_u -= 1
            if r_full_u < r_full_l:  # no full rules; [l,u] straddles a boundary
                result = 0.
            else:  # do the full rules
                n_l = self.starts[r_full_l]
                n_u = self.starts[r_full_u] + self.rules[r_full_u].npts - 1
                nodes = self.nodes[n_l:n_u+1]
                if self.node_dup[n_l] or self.node_dup[n_u]:
                    wts = self.wts[n_l:n_u+1].copy()
                    if self.node_dup[n_l]:
                        wts[0] = self.node_dup[n_l][1]
                    if self.node_dup[n_u]:
                        wts[-1] = self.node_dup[n_u][0]
                else:  # don't copy if we don't have to
                    wts = self.wts[n_l:n_u+1]
                # TODO:  If f() doesn't broadcast, iterate over nodes.
                result = sum(wts*f(nodes))

            # Add partial intervals at each end.
            if lower:
                range_method = getattr(self.rules[r_l], self.rm_name)
                result += range_method(f, lower[0], lower[1])
            if upper:
                range_method = getattr(self.rules[r_u], self.rm_name)
                result += range_method(f, upper[0], upper[1])

        return self.factor*result
