from sage.all import *  # type: ignore
from sage.graphs.distances_all_pairs import distances_all_pairs  # type: ignore
from sage.homology.chain_complex import ChainComplex  # type: ignore
from sage.matrix.constructor import matrix  # type: ignore
from sage.rings.all import ZZ  # type:ignore

BaseRing = ZZ


class MPSS:

    def __init__(self, graph, side, rmax=2, Reduced=True, Eulerian=False):

        self.graph = graph
        self.side = side
        self.rmax = rmax
        self.Reduced = Reduced
        self.Eulerian = Eulerian

        if self.rmax > 2:
            print("The specified rmax (3 or more) is not supported. The value will be adjusted to 2.")
            self.rmax = 2

        self.kmax = self.side + self.rmax
        self.lmax = int(self.side + self.rmax * (self.rmax - 1) / 2)

        self.dist = distances_all_pairs(self.graph)

        self.magnitude_complex = self._get_magnitude_complex()

        self._calculate_total_rank(0)

        if rmax == 0:
            return

        self.differentials = self._get_differentials()

        self.magnitude_homology = self._get_magnitude_homology()

        self._calculate_total_rank(1)

        if rmax == 1:
            return

        self.kmax -= 1

        self.integrated_magnitude_complex, self.integrated_magnitude_complex_reversed = self._get_integrated_magnitude_complex()
        self.basises = self._get_basises()
        self.images = self._get_images()
        self._shortening()
        self.differentials_dash = self._get_differentials_dash()
        self.magnitude_homology2 = self._get_magnitude_homology2()

        self.basises2 = self._get_basises2()
        self.images2 = self._get_images2()

        self._calculate_total_rank(2)

    def _get_magnitude_complex(self):

        MC_results = {(start, end, k, l): [] for start in self.graph.vertices() for end in self.graph.vertices() for k in range(self.kmax + 1) for l in range(self.lmax + 1)}

        if self.Eulerian:

            def add_generators(chain, l, adding_vertex, tracks, rests):
                k = len(chain) - 1
                if k <= self.kmax and l <= self.lmax:
                    MC_results[(chain[0], chain[-1], k, l)].append(chain)
                    for next_vertex in rests:
                        new_rests = rests.copy()
                        new_rests.remove(next_vertex)
                        new_tracks = tracks.copy()
                        new_tracks.add(next_vertex)
                        add_generators(chain + [next_vertex], l + self.dist[adding_vertex][next_vertex], next_vertex, new_tracks, new_rests)

            for initial_vertex in self.graph.vertices():
                all_vertices = set(self.graph.vertices())
                all_vertices.remove(initial_vertex)
                add_generators([initial_vertex], 0, initial_vertex, {initial_vertex}, all_vertices)

        else:

            def add_generators(chain, l, adding_vertex):
                k = len(chain) - 1
                if k <= self.kmax and l <= self.lmax:
                    MC_results[(chain[0], chain[-1], k, l)].append(chain)
                    for next_vertex in self.graph.vertices():
                        if next_vertex != adding_vertex:
                            add_generators(chain + [next_vertex], l + self.dist[adding_vertex][next_vertex], next_vertex)

            for initial_vertex in self.graph.vertices():
                add_generators([initial_vertex], 0, initial_vertex)

        for start in self.graph.vertices():
            for end in self.graph.vertices():
                for k in range(self.kmax + 1):
                    for l in range(self.lmax + 1):
                        MC_results[(start, end, k, l)] = {tuple(chain): i for (i, chain) in enumerate(MC_results[(start, end, k, l)])}

        return MC_results

    def _get_differentials(self):

        def differential(start, end, k, l):
            chains_here = self.magnitude_complex[(start, end, k, l)]
            another_chains = self.magnitude_complex[(start, end, k - 1, l)]
            boundary_matrix = {}
            for chain, i in chains_here.items():
                for z in range(len(chain) - 2):
                    if self.dist[chain[z]][chain[z + 1]] + self.dist[chain[z + 1]][chain[z + 2]] == self.dist[chain[z]][chain[z + 2]]:
                        j = another_chains[chain[: z + 1] + chain[z + 2 :]]
                        boundary_matrix[(j, i)] = boundary_matrix.get((j, i), 0) + (-1) ** (z + 1)
            return matrix(BaseRing, len(another_chains), len(chains_here), boundary_matrix)

        return {(start, end, k, l): differential(start, end, k, l) for start in self.graph.vertices() for end in self.graph.vertices() for k in range(1, self.kmax + 1) for l in range(self.lmax + 1)}

    def _get_magnitude_homology(self):

        def chains(start, end, l):
            differentials = {k: self.differentials[(start, end, k, l)] for k in range(1, self.kmax + 1) if self.magnitude_complex[(start, end, k, l)] or self.magnitude_complex[(start, end, k - 1, l)]}
            return ChainComplex(differentials, base_ring=BaseRing, degree=-1)

        def homology(start, end, l):
            return chains(start, end, l).homology(generators=True)

        return {(start, end, l): homology(start, end, l) for start in self.graph.vertices() for end in self.graph.vertices() for l in range(self.lmax + 1)}

    def _get_integrated_magnitude_complex(self):

        # Rank of magnitude complex at (k,l)
        self.total_length = {(k, l): 0 for k in range(self.kmax + 1) for l in range(self.lmax + 1)}
        # Dictionary of generators (in chain form), distinct from a basis
        generators2 = {(k, l): [] for k in range(self.kmax + 1) for l in range(self.lmax + 1)}
        # Reverse lookup dictionary
        generators2_rev = {(k, l): {} for k in range(self.kmax + 1) for l in range(self.lmax + 1)}
        # Helper for aggregating (start, end, k, l) results into (k, l) form
        self.segment = {(k, l): [0] for k in range(self.kmax + 1) for l in range(self.lmax + 1)}

        for k in range(self.kmax + 1):
            for l in range(self.lmax + 1):
                for start in self.graph.vertices():
                    for end in self.graph.vertices():
                        chains = self.magnitude_complex.get((start, end, k, l), [])
                        length = len(chains)
                        self.total_length[(k, l)] += length
                        self.segment[(k, l)].append(self.segment[(k, l)][-1] + length)
                        generators2[(k, l)].extend(chains)

                generators2_rev[(k, l)] = dict(enumerate(generators2[(k, l)]))
                generators2[(k, l)] = {tuple(chain): i for i, chain in generators2_rev[(k, l)].items()}

        return generators2, generators2_rev

    def _get_basises(self):

        basises = {(k, l): [] for k in range(self.kmax + 1) for l in range(self.lmax + 1)}
        # Aggregates (start, end, k, l) results into (k, l) form
        for k in range(self.kmax + 1):
            for l in range(self.lmax + 1):
                segment_phase = 0
                for start, end in ((s, e) for s in self.graph.vertices() for e in self.graph.vertices()):
                    homol = self.magnitude_homology.get((start, end, l))
                    if not homol:
                        segment_phase += 1
                        continue
                    for _, spaces_and_chains in homol.items():
                        if not spaces_and_chains:
                            continue
                        for space_and_chain in spaces_and_chains:
                            basis_vector = [0] * self.total_length[(k, l)]
                            basis = list(space_and_chain[1].vector(k))
                            place = self.segment[(k, l)][segment_phase]
                            for i, chain_coefficient in enumerate(basis):
                                basis_vector[place + i] = chain_coefficient
                            if any(basis_vector):
                                basises[(k, l)].append(basis_vector)
                    segment_phase += 1

        if self.Reduced:
            basises[(0, 0)] = []
            for i in range(len(self.integrated_magnitude_complex[(0, 0)]) - 1):
                v = [0] * len(self.integrated_magnitude_complex[(0, 0)])
                v[i] = -1
                v[i + 1] = 1
                basises[(0, 0)].append(v)

        return basises

    def _get_images(self):

        # Basis for the images of the differential matrices
        im_sep = {(start, end, k, l): [] for start in self.graph.vertices() for end in self.graph.vertices() for k in range(self.kmax + 1) for l in range(self.lmax + 1)}

        # Construct "im_sep"
        for k in range(self.kmax):
            for l in range(self.lmax + 1):
                for start in self.graph.vertices():
                    for end in self.graph.vertices():
                        M = list(self.differentials[(start, end, k + 1, l)])
                        if not any(M):
                            continue
                        rawl = len(M)
                        coll = len(M[0])
                        tM = [[0] * rawl for _ in range(coll)]
                        for i in range(coll):
                            for j in range(rawl):
                                tM[i][j] = M[j][i]
                        R = matrix(BaseRing, tM).rref()
                        for v in R:
                            if any(v):
                                im_sep[(start, end, k, l)].append(v)

        ims = {(k, l): [] for k in range(self.kmax + 1) for l in range(self.lmax + 1)}

        # Aggregates (start, end, k, l) results into (k, l) form
        for k in range(self.kmax + 1):
            for l in range(self.lmax + 1):
                segment_phase = 0
                for start, end in ((s, e) for s in self.graph.vertices() for e in self.graph.vertices()):
                    for im in im_sep[(start, end, k, l)]:
                        if im:
                            basis_vector = [0] * self.total_length[(k, l)]
                            basis = list(im)
                            place = self.segment[(k, l)][segment_phase]
                            for i, chain_coefficient in enumerate(basis):
                                basis_vector[place + i] = chain_coefficient
                            if any(basis_vector):
                                ims[(k, l)].append(basis_vector)
                    segment_phase += 1

        return ims

    def _shortening(self):
        # Subsequence of generators to be used for magnitude homology
        subseq_hom = {(k, l): set() for k in range(self.kmax + 1) for l in range(self.lmax + 1)}
        subseq_I = {(k, l): set() for k in range(self.kmax + 1) for l in range(self.lmax + 1)}
        subseq_K = {(k, l): set() for k in range(self.kmax + 1) for l in range(self.lmax + 1)}

        for k in range(self.kmax + 1):
            for l in range(self.lmax + 1):

                if not self.basises[(k, l)] and not self.images[(k, l)]:
                    continue

                vect_hom = [0] * len(self.integrated_magnitude_complex[(k, l)])
                for basis in self.basises[(k, l)]:
                    for i, coef in enumerate(basis):
                        if coef != 0:
                            vect_hom[i] = 1
                for i, v in enumerate(vect_hom):
                    if v != 0:
                        subseq_hom[(k, l)].add(i)

                vect_I = [0] * len(self.integrated_magnitude_complex[(k, l)])
                for im in self.images[(k, l)]:
                    for i, coef in enumerate(im):
                        if coef != 0:
                            vect_I[i] = 1
                for i, v in enumerate(vect_I):
                    if v != 0:
                        subseq_I[(k, l)].add(i)

                subseq_K[(k, l)] = subseq_hom[(k, l)] | subseq_I[(k, l)]

        # Rewriting with "subseq"
        for k in range(self.kmax + 1):
            for l in range(self.lmax + 1):
                if not subseq_K[(k, l)]:
                    self.integrated_magnitude_complex[(k, l)] = {}
                    self.integrated_magnitude_complex_reversed[(k, l)] = {}
                    continue
                newdict = {}
                for i, j in enumerate(subseq_K[(k, l)]):
                    newdict[tuple(self.integrated_magnitude_complex_reversed[(k, l)][j])] = i
                self.integrated_magnitude_complex[(k, l)] = newdict
                self.integrated_magnitude_complex_reversed[(k, l)] = {i: list(vect) for i, vect in enumerate(self.integrated_magnitude_complex[(k, l)])}

                for i, basis in enumerate(self.basises[(k, l)]):
                    newvect = [0] * len(self.integrated_magnitude_complex[(k, l)])
                    for j, aj in enumerate(subseq_K[(k, l)]):
                        newvect[j] = basis[aj]
                    self.basises[(k, l)][i] = newvect

                for i, im in enumerate(self.images[(k, l)]):
                    newvect = [0] * len(self.integrated_magnitude_complex[(k, l)])
                    for j, aj in enumerate(subseq_K[(k, l)]):
                        newvect[j] = im[aj]
                    self.images[(k, l)][i] = newvect

    def _get_differentials_dash(self):

        def differential_dash(k, l):
            chains_here = self.basises.get((k, l), [])
            another_chains = self.basises.get((k - 1, l - 1), []) + self.images.get((k - 1, l - 1), [])
            sep = len(self.basises.get((k - 1, l - 1), []))
            boundary_matrix = []
            for basis in chains_here:
                destination = [0] * len(self.integrated_magnitude_complex[(k - 1, l - 1)])
                for chain_number, chain_coefficient in enumerate(basis):
                    if chain_coefficient == 0:
                        continue
                    chain = self.integrated_magnitude_complex_reversed[(k, l)][chain_number]
                    if self.dist[chain[0]][chain[1]] == 1:
                        j = self.integrated_magnitude_complex[(k - 1, l - 1)][tuple(chain[1:])]
                        destination[j] += chain_coefficient
                    for z in range(len(chain) - 2):
                        if self.dist[chain[z]][chain[z + 1]] + self.dist[chain[z + 1]][chain[z + 2]] == self.dist[chain[z]][chain[z + 2]] + 1:
                            j = self.integrated_magnitude_complex[(k - 1, l - 1)][tuple(chain[: z + 1] + chain[z + 2 :])]
                            destination[j] += (-1) ** (z + 1) * chain_coefficient
                    if self.dist[chain[-2]][chain[-1]] == 1:
                        j = self.integrated_magnitude_complex[(k - 1, l - 1)][tuple(chain[:-1])]
                        destination[j] += (-1) ** (len(chain) - 1) * chain_coefficient
                # represent "destination" in the basis of "another_chains"
                basis_matrix = matrix(BaseRing, len(another_chains), len(destination), another_chains)
                coefficients = list(basis_matrix.solve_left(matrix(BaseRing, 1, len(destination), destination))[0])
                boundary_matrix.append(coefficients[:sep])
            return matrix(BaseRing, len(chains_here), len(self.basises.get((k - 1, l - 1), [])), boundary_matrix).transpose()

        return {(k, l): differential_dash(k, l) for k in range(1, self.kmax + 1) for l in range(k, k + self.lmax + 1) if self.basises.get((k, l)) or self.basises.get((k - 1, l - 1))}

    def _get_magnitude_homology2(self):

        def chains_dash(l_minus_k):
            differentials = {k: self.differentials_dash[(k, l_minus_k + k)] for k in range(1, self.kmax + 1) if self.basises.get((k, l_minus_k + k)) or self.basises.get((k - 1, l_minus_k + k - 1))}
            return ChainComplex(differentials, base_ring=BaseRing, degree=-1)

        def homology_dash(l_minus_k):
            return chains_dash(l_minus_k).homology(generators=True)

        return {l_minus_k: homology_dash(l_minus_k) for l_minus_k in range(self.lmax + 1)}

    def _get_basises2(self):

        # Homology generators, represented as coefficient vectors for generators2
        basises = {(k, l): [] for k in range(self.kmax + 1) for l in range(self.lmax + 1)}

        # Construct "basises"
        for k in range(self.kmax + 1):
            for l in range(self.lmax + 1):
                homol = self.magnitude_homology2.get((l - k))
                if not homol:
                    continue
                for _, spaces_and_chains in homol.items():
                    if not any(spaces_and_chains):
                        continue
                    for space_and_chain in spaces_and_chains:
                        chain = list(space_and_chain[1].vector(k))
                        v = [0] * len(self.integrated_magnitude_complex[(k, l)])
                        for i, chain_coefficient in enumerate(chain):
                            if chain_coefficient == 0:
                                continue
                            add_vect = self.basises[(k, l)][i]
                            for i in range(len(add_vect)):
                                v[i] += chain_coefficient * add_vect[i]
                        if any(v):
                            basises[(k, l)].append(v)

        return basises

    def _get_images2(self):

        ims = {(k, l): [] for k in range(self.kmax + 1) for l in range(self.lmax + 1)}

        # Construct "ims"
        for k in range(self.kmax):
            for l in range(self.lmax - 1):
                M = list(self.differentials_dash.get((k + 1, l + 1), []))
                if any(M):
                    rawl = len(M)
                    coll = len(M[0])
                    tM = [[0] * rawl for _ in range(coll)]
                    for i in range(coll):
                        for j in range(rawl):
                            tM[i][j] = M[j][i]
                    R = matrix(BaseRing, tM).rref()
                    for chain in R:
                        v = [0] * len(self.integrated_magnitude_complex[(k, l)])
                        for i, chain_coefficient in enumerate(chain):
                            if chain_coefficient == 0:
                                continue
                            add_vect = self.basises[(k, l)][i]
                            for i in range(len(add_vect)):
                                v[i] += chain_coefficient * add_vect[i]
                        if any(v):
                            ims[(k, l)].append(v)

        return ims

    def _calculate_total_rank(self, r):

        def myrank(I):
            c = 0
            for v in I:
                if v == 0:
                    c += 1
            return c

        if r == 0:
            mcpx = self.magnitude_complex
            self.total_rank = [{(k, l): 0 for k in range(self.side + 1) for l in range(self.side + 1)}]

            for start, end, k, l in ((start, end, k, l) for start in self.graph.vertices() for end in self.graph.vertices() for k in range(self.side + 1) for l in range(self.side + 1)):
                self.total_rank[0][(k, l)] += len(mcpx[(start, end, k, l)])

        elif r == 1:
            if self.rmax < 1:
                return
            homology = self.magnitude_homology
            self.total_rank.append({(k, l): 0 for k in range(self.side + 1) for l in range(self.side + 1)})

            for start, end, l, degree, group in (
                (start, end, l, degree, group_and_chain[0])
                for start in self.graph.vertices()
                for end in self.graph.vertices()
                for l in range(self.side + 1)
                for degree, group_and_chains in sorted(homology[(start, end, l)].items())
                for group_and_chain in group_and_chains
            ):
                if degree < 0:
                    continue
                self.total_rank[1][degree, l] += myrank(group.invariants())

            if self.Reduced:
                self.total_rank[1][(0, 0)] -= 1

        elif r == 2:
            if self.rmax < 2:
                return
            homology2 = self.magnitude_homology2

            self.total_rank.append({(k, l): 0 for k in range(self.side + 1) for l in range(self.side + 1)})

            for l_minus_k, degree, group in (
                (l_minus_k, degree, group_and_chain[0])
                for l_minus_k in range(self.side + 1)
                for degree, group_and_chains in sorted(homology2[l_minus_k].items())
                for group_and_chain in group_and_chains
            ):
                if degree > self.side or degree < 0 or l_minus_k + degree > self.side:
                    continue
                self.total_rank[2][degree, l_minus_k + degree] += myrank(group.invariants())

    #################################################################################################################################################

    def print_sheet(self, r, twidth=6):

        if r == 0:

            print("0th sheet (Magnitude Complex):\n l\\k", end="")
            print("".join(f"{l:{twidth}d}" for l in range(self.side + 1)))

            for l in range(self.side + 1):
                print(f"{l:4d}", end="")
                for k in range(self.side + 1):
                    value = self.total_rank[0][k, l]
                    print(f"{value:{twidth}d}" if value != 0 else f"{'':{twidth-1}}-", end="")
                print()

        elif r == 1:

            if self.rmax < 1:
                print("Error: r is too large")
                return

            print("1st sheet (Magnitude Homology):\n l\\k", end="")
            print("".join(f"{l:{twidth}d}" for l in range(self.side + 1)))

            for l in range(self.side + 1):
                print(f"{l:4d}", end="")
                for k in range(self.side + 1):
                    value = self.total_rank[1][(k, l)]
                    print(f"{value:{twidth}d}" if value != 0 else f"{'':{twidth-1}}-", end="")
                print()

        elif r == 2:

            if self.rmax < 2:
                print("Error: r is too large")
                return

            print("2nd sheet(Bigraded Path Homology):\n l\\k", end="")
            print("".join(f"{l:{twidth}d}" for l in range(self.side + 1)))

            for l in range(self.side + 1):
                print(f"{l:4d}", end="")
                for k in range(self.side + 1):
                    value = self.total_rank[2][(k, l)]
                    print(f"{value:{twidth}d}" if value != 0 else f"{'':{twidth-1}}-", end="")
                print()

    def get_generators(self, k, l, r):

        if max(k, l) > self.side:
            print("Error: k or l is too large")
            return

        def chain_form(k, l, vect, r=1):
            result = ""
            for chain_number, chain_coefficient in enumerate(vect):
                if chain_coefficient != 0:
                    chain = self.integrated_magnitude_complex_reversed[(k, l)][chain_number]
                    if chain_coefficient != 1:
                        result += f"{chain_coefficient}"
                    result += f"{chain}+"
            result = result[:-1]
            result = result.replace("+-", "-")
            result = result.replace("-1", "-")
            return result

        results = []

        if r == 0:
            for start in self.graph.vertices():
                for end in self.graph.vertices():
                    for chain in self.magnitude_complex[(start, end, k, l)]:
                        results.append(str(list(chain)))

        elif r == 1:
            if self.rmax < 1:
                print("Error: r is too large")
                return
            for vect in self.basises[(k, l)]:
                if chain_form(k, l, vect) != "":
                    results.append(chain_form(k, l, vect))

        elif r == 2:
            if self.rmax < 2:
                print("Error: r is too large")
                return
            for vect in self.basises2[(k, l)]:
                if chain_form(k, l, vect) != "":
                    results.append(chain_form(k, l, vect))

        return results

    def get_divisors(self, k, l, r=1):

        if max(k, l) > self.side:
            print("Error: k or l is too large")
            return

        def chain_form(k, l, vect, r=1):
            result = ""
            for chain_number, chain_coefficient in enumerate(vect):
                if chain_coefficient != 0:
                    chain = self.integrated_magnitude_complex_reversed[(k, l)][chain_number]
                    if chain_coefficient != 1:
                        result += f"{chain_coefficient}"
                    result += f"{chain}+"
            result = result[:-1]
            result = result.replace("+-", "-")
            result = result.replace("-1", "-")
            return result

        results = []

        if r == 0:
            return results

        elif r == 1:
            if self.rmax < 1:
                print("Error: r is too large")
                return
            for vect in self.images[(k, l)]:
                if chain_form(k, l, vect) != "":
                    results.append(chain_form(k, l, vect))

        elif r == 2:
            if self.rmax < 2:
                print("Error: r is too large")
                return
            for vect in self.images2[(k, l)]:
                if chain_form(k, l, vect) != "":
                    results.append(chain_form(k, l, vect))

        return results

    # def get_differential(self, k, l, r, start=None, end=None):
    #    return

    def show_space_info(self, k, l, r):

        if max(k, l) > self.side:
            print("Error: k or l is too large")
            return

        if r == 0:
            result_length = 0
            for start in self.graph.vertices():
                for end in self.graph.vertices():
                    cpx = self.magnitude_complex.get((start, end, k, l))
                    result_length += len(cpx)
            print(f"({k}, {l}, {r}): Z^{result_length}")
        elif r == 1:
            if self.rmax < 1:
                print("Error: r is too large")
                return
            result = ""
            for start in self.graph.vertices():
                for end in self.graph.vertices():
                    homol = self.magnitude_homology.get((start, end, l))
                    if not homol:
                        continue
                    spaces_and_chains = homol.get(k)
                    if not spaces_and_chains:
                        continue
                    for space_and_chain in spaces_and_chains:
                        result += str(space_and_chain[0]) + "×"
            if not result:
                print(f"({k}, {l}, {r}): {0}")
                return
            result = result[:-1]
            if self.Reduced and (k, l) == (0, 0):
                result = result[:-2]
            print(f"({k}, {l}, {r}): {result}")

        elif r == 2:
            if self.rmax < 2:
                print("Error: r is too large")
                return
            result = ""
            homol = self.magnitude_homology2.get(l - k)
            if not homol:
                print(f"({k}, {l}, {r}): {0}")
                return
            spaces_and_chains = homol.get(k)
            if not spaces_and_chains:
                print(f"({k}, {l}, {r}): {0}")
                return
            for space_and_chain in spaces_and_chains:
                result += str(space_and_chain[0]) + "×"
            if not result:
                print(f"({k}, {l}, {r}): {0}")
                return
            result = result[:-1]
            if self.Reduced and (k, l) == (0, 0):
                result = result[:-2]
            print(f"({k}, {l}, {r}): {result}")

    def get_homology_group(self, l, r, start=None, end=None):  # 未完成

        if l > self.side:
            print("Error: l is too large")
            return

        if r == 0:
            if not start or not end:
                print("Error: start and end must be specified")
                return
            return {k: self.magnitude_complex.get((start, end, k, l)) for k in range(self.side + 1)}

        if r == 1:
            if self.rmax < 1:
                print("Error: r is too large")
                return
            if not start or not end:
                print("Error: start and end must be specified")
                return
            return self.magnitude_homology.get((start, end, l))

        if r == 2:
            if self.rmax < 2:
                print("Error: r is too large")
                return
            result = {k: [] for k in range(self.side + 1)}
            for k in range(self.side + 1):
                result[k] = self.magnitude_homology2.get(l - k)
            return result

    def magnitude(self):

        if self.rmax < 1:
            print("Error: r is too large")
            return

        mag = ""

        for l in range(self.side + 1):
            coef = 0
            for k in range(self.side + 1):
                coef += (-1) ** k * self.total_rank[1][(k, l)]
            if coef != 0:
                if coef > 0 and mag != "":
                    mag += "+"
                if coef == -1:
                    mag += "-"
                elif coef != 1:
                    mag += str(coef)
                if l == 0:
                    mag += str(l)
                else:
                    mag += f"q^{l}"

        if mag != "":
            mag += "+⋯⋯"

        return mag
