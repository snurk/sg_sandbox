#!/usr/bin/env python3
import sys
import gfapy

if len(sys.argv) < 3:
  print("Usage: %s <gfa in> <gfa out> [name prefix (default: 'm_')]" % sys.argv[0])
  sys.exit(239)

print("Loading graph")
gfa = gfapy.Gfa.from_file(sys.argv[1], version = "gfa1")
print("Loaded")

name_prefix = "m_"
if len(sys.argv) > 3:
  name_prefix = sys.argv[3]
#for p in gfa.linear_paths():
#  print("Will merge path: ", p)

def create_merged_segment(gfa, segpath,
    merged_name=None, enable_tracking=False, cut_counts=False):
  merged = gfa.try_get_segment(segpath[0].segment).clone()
  merged_vlevel = merged.vlevel
  merged.vlevel = 0
  total_cut = 0
  a = segpath[0]
  first_reversed = (a.end_type == "L")
  last_reversed = None
  if merged_name == "short":
    merged_name = gfa.unused_name()
  gfa._add_segment_to_merged(merged, gfa.segment(a.segment),
      first_reversed, 0, True, enable_tracking=enable_tracking,
      merged_name=merged_name)
  #for i in range(len(segpath)-1):
  #  b = gfapy.SegmentEnd(segpath[i+1]).inverted()
  for s in segpath[1:]:
    b = gfapy.SegmentEnd(s).inverted()
    ls = gfa.segment(a.segment).end_relations(a.end_type, b, "dovetails")
    if len(ls) != 1:
      msg = "A single link was expected between {}".format(a) + \
            "and {}".format(b) + "{} were found".format(len(ls))
      raise gfapy.ValueError(msg)
    l = ls[0]
    if not l.overlap:
      cut = 0
    else:
      cut = min(l.overlap.length_on_query(), gfa.segment(b.segment).LN)
    #elif all(op.code in ["M","="] for op in l.overlap):
    #  cut = sum([len(op) for op in l.overlap])
    #else:
    #  raise gfapy.ValueError(
    #      "Merging is only allowed if all operations are M/=")
    total_cut += cut
    last_reversed = (b.end_type == "R")
    gfa._add_segment_to_merged(merged, gfa.segment(b.segment),
        last_reversed, cut, False, enable_tracking=enable_tracking,
        merged_name=merged_name)
    a = gfapy.SegmentEnd(b).inverted()
  merged.vlevel = merged_vlevel
  if isinstance(merged.name, list):
    merged.name = "_".join(merged.name)
  ortag = merged.get("or")
  if isinstance(ortag, list):
    merged.set("or", ",".join(ortag))
  if not gfapy.is_placeholder(merged.sequence):
    merged.sequence = "".join(merged.sequence)
    if not merged.LN:
      merged.LN = len(merged.sequence)
    elif gfa._vlevel > 0 and merged.LN != len(merged.sequence):
      raise gfapy.InconsistencyError(
          "Computed sequence length {} ".format(merged.sequence.length)+
          "and computed LN {} differ".format(merged.LN))
  if merged.length is not None:
    for count_tag in ["KC", "RC", "FC"]:
      merged.set(count_tag, None)
  else:
    factor = 1
    if cut_counts:
      factor = merged.length / (total_cut+merged.length)
    for count_tag,count in gfa.__sum_of_counts(segpath,factor).items():
      merged.set(count_tag, count)
  return merged, first_reversed, last_reversed

def link_merged(gfa, merged_name, segment_end, is_reversed):
  #print("")
  #print("Linking to", segment_end, " ", is_reversed)
  #print("Details", segment_end.segment, " ", segment_end.end_type)
  for l in gfa.segment(segment_end.segment).dovetails_of_end(segment_end.end_type).copy():
    #print("Linking ", l)
    l2 = l.clone()
    if l2.to_segment == segment_end.segment:
      l2.to_segment = merged_name
      if is_reversed:
        l2.to_orient = gfapy.invert(l2.to_orient)
    else:
      l2.from_segment = merged_name
      if is_reversed:
        l2.from_orient = gfapy.invert(l2.from_orient)
    l.disconnect()
    gfa.add_line(l2)

def merge_linear_path(gfa, segpath,
                      enable_tracking=False, merged_name=None,
                      cut_counts=False):
  """Merge a specified linear path of dovetail overlaps connecting segments.
  Note:
    for the parameter usage, see merge_linear_paths();
    the only difference is that merged_name can be set to a string (different
    from 'short'), which will be used as a name for the merged segment.
  """
  if len(segpath) < 2:
    return gfa
  segpath = [gfapy.SegmentEnd(s) for s in segpath]
  #print("Merging path", segpath)
  merged, first_reversed, last_reversed = create_merged_segment(gfa, segpath,
          merged_name=merged_name,cut_counts=cut_counts,
          enable_tracking=enable_tracking)
  gfa.append(merged)
  link_merged(gfa, merged.name, segpath[0].inverted(), first_reversed)
  link_merged(gfa, merged.name, segpath[-1], last_reversed)
  idx1 = 0
  idx2 = None
  for sn_et in segpath[idx1:idx2]:
    gfa.segment(sn_et.segment).disconnect()
  #return gfa

  """Find and merge linear paths of dovetail overlaps connecting segments.
  Note:
    Besides obviously the dovetail overlaps, all lines refererring to the
    merged segments (containments, internal edges, paths, sets, fragments,
    gaps) are removed from the Gfa instance.
  Parameters:
    merged_name (str): if 'short', then a name is computed using an unused
      integer; otherwise the name is computed using a combination of the
      names of the merged segments, separated by an underscore
    cut_counts (bool): if True, the total count in merged segment m,
       composed of segments s of set S is multiplied by the factor
       ``Sum(|s in S|)/|m|``
    enable_tracking: if True, tracking information is added as follows;
      the name of the component segments is stored in the ``or`` tag (or the
      content of their ``or`` tag, instead of the name, if any) and their
      starting positions is stored in the ``mp`` tag; the ``rn`` tag, used
      for storing possibe inversion positions by the random orientation
      methods of this library, is inherited and the positions updated;
      unless merged_name is set to 'short', the computation of the merged
      name is enhanced, in that reverse complement components are suffixed
      with ``^`` and parenthesis added by the random orientation methods of
      this library are inherited
  """

def end_type_direction_char(end_type):
  if end_type == "L":
    return "-"
  elif end_type == "R":
    return "+"
  assert(False)

def merge_linear_paths(gfa,
                      name_prefix,
                      enable_tracking=False,
                      cut_counts=False):
  print("Getting linear paths")
  paths = gfa.linear_paths()
  print("Collected")
  psize = sum([len(path) for path in paths])
  print("Merging %d paths (total %d segments)" % (len(paths), psize))
  for i, path in enumerate(paths):
    #print("Processing unipath #%d consisting of %d nodes" % (i, len(path)))
    name = "%s%d" % (name_prefix, i)
    merge_linear_path(gfa, path,
                      merged_name=name,
                      cut_counts=cut_counts,
                      enable_tracking=enable_tracking)
    print(name, ','.join([str(s.segment) + end_type_direction_char(s.end_type) for s in path]), file=sys.stderr)
  #return gfa

merge_linear_paths(gfa, name_prefix, enable_tracking=False)

print("Storing to file", sys.argv[2])

gfa.to_file(sys.argv[2])

print("Compression done")
