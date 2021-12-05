from patch import p
import arbor
import itertools

class VirtualSection:
    def __init__(self, stitches):
        pass


def _cocv(data, attr):
    yield from enumerate(map(getattr(data, attr), range(data.num_cv())))


def guess_subcv_precision(cc, label, cv):
    region = f'"{label}"'
    cv_data = arbor.cv_data(cc)
    cv_cables = cv_data.cables(cv)
    region_as_cables = cc.cables(region)
    region_as_cvs = arbor.intersect_region(cc, region, cv_data)
    parent = cv_data.parent(cv)
    print("CV cables:", cv, cv_cables)
    print("Problematic region cables:", label, region_as_cables)
    print("Problematic region cv proportions:", label, region_as_cvs)
    print("Gogo, solve it Robin")
    if len(region_as_cables) == 1:
        if len(region_as_cvs) == 1:
            print("Decided to solve it by claiming only cable", region_as_cables[0])
            return region_as_cables[0]
        else:
            subs = [i for i, r in enumerate(region_as_cvs) if r[1] != 1]
            print("incomplete subs", subs)
            if len(subs) == 1:
                print("Decided to solve it by associating only noncomplete cv", region_as_cvs[subs[0]], "to only cable", region_as_cables[0])
                print("Returning cable", region_as_cables[0])
                return region_as_cables[0]
            raise RuntimeError("SubCV solver expected region cables to match region CVs")
        return region_as_cables[0]
    raise RuntimeError(f"Couldn't solve subcv arbor instruction")


def arb_to_sections(morphology, label_dict, decor):
    print(dir(morphology))
    cc = arbor.cable_cell(morphology, label_dict, decor)
    cv_data = arbor.cv_data(cc)
    # for i, cable in _cocv(cv_data, "cables"):
    #     print(i, cable)
    cv_labels = {}
    cv_splits = {}
    cv_split_labels = {}
    for label in label_dict:
        if label == "all" or "(" in label:
            # Skip `all` label and anti-overpainting constructs.
            continue
        for cv, prop in arbor.intersect_region(cc, f'"{label}"', cv_data):
            if prop == 1:
                cv_labels.setdefault(cv, set()).add(label)
            else:
                cable = guess_subcv_precision(cc, label, cv)
                cv_splits.setdefault(cv, []).append(cable)
                cv_split_labels.setdefault(cv, []).append(label)
        # print(label, list((cv_cables[a[0]], a[1]) for a in r))
    # print(dir(cv_data))
    # print(cv_labels)
    # print(cv_splits)
    tree, cables, skipped_parents = collect_cv_tree(cv_data)
    print(tree)
    print(cables)
    new_id = itertools.count(cv_data.num_cv())
    replaced_ends = {}
    for cv, splits, labels in zip(cv_splits.keys(), cv_splits.values(), cv_split_labels.values()):
        split_cv(tree, cables, cv_data, cv, splits, cv_labels, labels, replaced_ends, new_id, skipped_parents)
    print(tree)
    print(cables)
    # Create sections from the tree
    ptr = tree[-1][0]
    print("starting on", ptr)
    sec = []
    sec_id = 0
    seclist = []
    sec_tree = {sec_id: []}
    ptr_stack = []
    parent_stack = []
    parent = None

    def finish_sec(nxt):
        nonlocal ptr, sec, sec_id, seclist, sec_tree, parent
        seclist.append(sec)
        if parent is not None:
            sec_tree.setdefault(parent, []).append(sec_id)
        parent = sec_id
        sec_id += 1
        sec = []
        ptr = nxt

    print(*_cocv(cv_data, "cables"))
    print([len(c[1]) for c in _cocv(cv_data, "cables")])
    while True:
        sec.append(ptr)
        children = tree.get(ptr, [])
        print("going for children", children)
        if not children:
            try:
                ptr = ptr_stack.pop()
                finish_sec(ptr)
                if not parent_stack:
                    parent = None
                else:
                    parent = parent_stack.pop()
            except IndexError:
                break
        elif len(children) == 1:
            if cv_labels[ptr] != cv_labels[children[0]]:
                finish_sec(children[0])
            else:
                ptr = children[0]
        else:
            print("Multiple children!")
            finish_sec(children[0])
            ptr_stack.extend(children[1:])
            parent_stack.extend([sec_id] * len(children[1:]))
    print(seclist)
    print([[cv_labels[cv] for cv in cvs] for cvs in seclist])



def split_cv(tree, cables, cv_data, cv, splits, labels, split_labels, replaced_ends, new_id, skipped_parents):
    unsplit_labels = labels[cv]
    parts = list(zip(new_id, sorted(splits, key=lambda x: x.prox)))
    start = parts[0][0]
    # Claim cables for the new CVs
    for new_cv, cable in parts:
        cables[new_cv] = cable
        labels[new_cv] = unsplit_labels.copy()
    for (new_cv, _), split_label in zip(parts, split_labels):
        labels[new_cv].add(split_label)
    # Connect new subsequent pieces to the next
    for (a, ca), (b, cb) in pairwise(parts):
        if ca.dist != cb.prox:
            raise NotImplementedError("Sorry, still have to deal with gaps, open an issue.")
        tree[a] = [b]
    # Connect the children of the unsplit CV to the last new CV
    tree[b] = tree[cv]
    del tree[cv]
    replaced_ends[cv] = b
    # Connect the first child to the parent of the unsplit CV
    parent = _skip_parents(cv_data.parent(cv), skipped_parents)
    true_parent = replaced_ends.get(parent, parent)
    parent_children = tree[true_parent]
    parent_children.remove(cv)
    parent_children.append(start)

def collect_cv_tree(cv_data):
    tree = {}
    skipped_parents = {}
    for cv, parent in _cocv(cv_data, "parent"):
        if len(cv_data.cables(cv)) > 1:
            skipped_parents[cv] = parent
        else:
            while parent in skipped_parents:
                parent = skipped_parents[parent]
            tree.setdefault(parent, []).append(cv)
    cables = {cv: cable for cv, cable in _cocv(cv_data, "cables")}
    return tree, cables, skipped_parents


def _skip_parents(parent, skippers):
    while parent in skippers:
        parent = skippers[parent]
    return parent


def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def _dims(segments):
    # Get the (prox - dist) ^ 2 along a certain axis
    d = lambda seg, axis: (getattr(seg.prox, axis) - getattr(seg.dist, axis)) ** 2
    # Sum and sqrt the distances along the x, y, z axes for eucl. dist
    eucl = [sum(d(s, axis) for axis in ("x", "y", "z")) ** (1/2) for s in segments]
    # Average the frustrum prox and dist radius for equivalent segment cilinder radii
    radii = [(s.prox.radius + s.dist.radius) / 2 for s in segments]
    total_length = sum(eucl)
    # Average of the segment cilinder radii, weighted by segment length
    r = sum(r * seg_len for r, seg_len in zip(radii, eucl)) / total_length / len(radii)
    # Return `L` and `diam`
    return total_length, r * 2


def make_section(segments):
    import neuron
    sec = neuron.h.Section()
    sec.L, sec.diam = _dims(segments)
    return sec

def make_cell(morphology, label_dict, decor):
    cell = NrnCell()
    num = morphology.num_branches
    cell.all = [
        make_section(morphology.branch_segments(i))
        for i in range(num)
    ]
    for i, sec in enumerate(cell.all):
        p = morphology.branch_parent(i)
        # If branch has no parent, this beaut is returned
        # from the `branch_parent` function
        if p != 4294967295:
            sec.connect(cell.all[p])

    return cell
