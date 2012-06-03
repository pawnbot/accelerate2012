#include <string>
#include <algorithm>
#include <map>
#include <iostream>
#include <fstream>

#include "solver.h"
#include "prefix_tree.h"

using namespace::std;

typedef std::pair<int, int> suffix_pair;

#define MAX_PART_LEN 10

static void parse_file (
  const char* filename,
  string_system* my_system,
  bool referal
  )
{
  ifstream stream (filename,ios::in);
  string name = "";
  string dna_sequence = "";

  while (!stream.eof ())
    {
      string new_str;
      getline (stream, new_str);
      if (new_str[0] == '>')
        {
          if (!name.empty ())
            {
              if (!referal)
                my_system->add_string (name, dna_sequence);
              else
                my_system->set_referal (dna_sequence);
              dna_sequence.clear ();
            }

          name = new_str.substr (1);
        }
      else
        {
          dna_sequence += new_str;
        }
    }
  if (!name.empty ())
    {
      if (!referal)
        my_system->add_string (name, dna_sequence);
      else
        my_system->set_referal (dna_sequence);
    }
}

static void find_matches (int thread_index, int thread_count, int moves_count, prefix_tree* tree,
                          const string &ref, const string &str, bool swap_pos, vector<matches_pair> &matches)
{
  int part_len = min (moves_count, MAX_PART_LEN);
  int shift = 1;
  for (int i = 1; i < part_len; ++i)
    shift = (shift<<2);

  size_t str_length = str.size ();

  map<int, int> *old_layer = new map<int, int> ();
  map<int, int> *new_layer = new map<int, int> ();
  map<int, int>::iterator it;

  matches_pair new_match;

  int first = thread_index * str_length / thread_count;
  int last = (thread_index + 1) * str_length / thread_count;
  if ((size_t) last == str_length)
    last--;
  int start_pos = max (first, moves_count - 1);
  int part_str = 0;
  int node, pos, len;
  int buf, prev_node = -1, prev_len = -INFINITY;
  int depth_in_edge = 0;
  char prev_char = '\0';
  int prev_i = last;


  for (int i = last; i > max (0, last - part_len + 1); --i)
    part_str = (part_str>>2) + (char_code (str[i]) - 1) * shift;

  string::const_iterator next_char = str.begin () + max (0, last - part_len + 1);
  for (int i = last; i >= start_pos; --i, --next_char)
    {
      part_str = (part_str>>2) + (char_code (*next_char) - 1) * shift;
      if (tree->part_index (part_str) < 0)
        {
          continue;
        }
      else
        {
          if (prev_i - i > 1)
            {
              prev_len = -INFINITY;
              old_layer->clear ();
            }
          if (prev_len + depth_in_edge <= part_len)
            {
              node = tree->part_index (part_str);
              pos  = tree->part_pos (part_str);
              len = part_len;
            }
          else
            {
              prev_len--;
              prev_node = tree->link_jump (prev_node);
              len = prev_len;
              tree->edge_jump (prev_char, prev_node, depth_in_edge, node, pos, len);
            }

          string::const_iterator str_char = str.begin () + i - len;
          for (; len < moves_count; ++len, --str_char)
            {
              buf = node;
              tree->move (*str_char, node, pos);
              if (node < 0)
                break;
              if (buf != node)
                {
                  prev_node = buf;
                  prev_len = len;
                  depth_in_edge = 0;
                  prev_char = *str_char;
                }
              depth_in_edge++;
            }

          if (node >= 0)
            {
              vector<int> indexes;
              tree->get_indexes (node, indexes);
              for (size_t k = 0; k < indexes.size (); ++k)
                {
                  int j = indexes[k];
                  int result = len;
                  if (old_layer->count (j + 1))
                    result = old_layer->at (j + 1) - 1;
                  else
                    {
                      while (i - result >= 0 && j - result >= 0 &&
                             str[i - result] == ref[j - result])
                        result++;
                    }
                  new_layer->insert (suffix_pair (j, result));
                }
              if (i < last || i + 1 == (int) str_length)
                for (it = new_layer->begin (); it != new_layer->end (); ++it)
                  {
                    int j = it->first;
                    int result = it->second;
                    if (old_layer->count (j + 1))
                      continue;
                    if (old_layer->count (j) && result <= old_layer->at (j))
                      continue;
                    if (new_layer->count (j + 1) && result <= new_layer->at (j + 1))
                      continue;
                    if (swap_pos)
                      {
                        new_match.ref_start = j - result + 2;
                        new_match.ref_end = j + 1;
                        new_match.str_start = i - result + 2;
                        new_match.str_end = i + 1;
                      }
                    else
                      {
                        new_match.ref_start = i - result + 2;
                        new_match.ref_end = i + 1;
                        new_match.str_start = j - result + 2;
                        new_match.str_end = j + 1;
                      }
                    matches.push_back (new_match);
                  }
              prev_i = i;
              swap (old_layer, new_layer);
              new_layer->clear ();
            }
        }
    }

  FREE (old_layer);
  FREE (new_layer);
}

void solve (
  size_t thread_index,
  size_t thread_count,
  int argc,
  char* argv[],
  string_system* my_system,
  pthread_barrier_t *sync
  )
{
  int min_match_length = atol (argv[2]);
  int part_len = min (min_match_length, MAX_PART_LEN);

  for (int i = 3; i < argc; ++i)
    {
      if (i % thread_count == thread_index)
        {
          if (i == 3)
            parse_file (argv[i], my_system, true);
          else
            parse_file (argv[i], my_system, false);
        }
    }
  pthread_barrier_wait (sync);
  /// Reading done

  size_t length = my_system->ref.size ();
  const string &ref = my_system->ref;

  static size_t total_len = 0;

  if (thread_index == 0)
    {
      /// TODO: Make parallel
      sort (my_system->data.begin (), my_system->data.end ());
      for (size_t i = 0; i < my_system->data.size (); ++i)
        total_len += my_system->data[i].str.size ();
    }
  pthread_barrier_wait (sync);

  static prefix_tree *ref_tree = 0;
  static prefix_tree *str_tree = 0;

  if (length < total_len && thread_index == 0)
    {
      ref_tree = new prefix_tree (part_len);
      ref_tree->parse (ref);
    }

  pthread_barrier_wait (sync);

  if (ref_tree)
    {
      ref_tree->calculate_depth (thread_index, thread_count, sync);
    }

  pthread_barrier_wait (sync);
  /// Referal string now in memory
  size_t size = my_system->data.size ();

  for (size_t index = thread_index; index < size; index += thread_count)
    {
      const string &str = my_system->data[index].str;
      size_t str_length = str.size ();

      if (length > 10 && str_length > 10 && length * str_length > 5000000)
        continue;
      /// Applying this alg only on tiny strings
      size_t i, j;
      char new_char;
      matches_pair new_match;
      vector<matches_pair> matches;
      int* old_layer = new int [str_length + 1];
      int* new_layer = new int [str_length + 1];
      new_char = ref[0];
      for (j = 0; j < str_length; ++j)
        old_layer[j] = (str[j] == new_char);
      for (i = 1; i < length; ++i)
        {
          new_char = ref[i];
          new_layer[0] = (str[0] == new_char);
          for (j = 1; j < str_length; ++j)
            new_layer[j] = (old_layer[j - 1] + 1) * (str[j] == new_char);

          for (j = 0; j < str_length; ++j)
            {
              if (old_layer[j] >= min_match_length)
                {
                  //this match can be shifted on i and j
                  if (j + 1 < str_length && old_layer[j] <= new_layer[j + 1])
                    continue;
                  //this match can be shifted on j
                  if (j + 1 < str_length && old_layer[j] <= old_layer[j + 1])
                    continue;
                  //this match can be shifted on i
                  if (old_layer[j] <= new_layer[j])
                    continue;
                  new_match.ref_start = i - old_layer[j] + 1;
                  new_match.ref_end = i;
                  new_match.str_start = j - old_layer[j] + 2;
                  new_match.str_end = j + 1;
                  matches.push_back (new_match);
                }
            }
          swap (old_layer, new_layer);
        }
      i = length;
      for (j = 0; j < str_length; ++j)
        {
          if (old_layer[j] >= min_match_length)
            {
              //this match can be shifted on i and j
              /// FALSE
              //this match can be shifted on j
              if (j + 1 < str_length && old_layer[j] <= old_layer[j + 1])
                continue;
              //this match can be shifted on i
              /// FALSE
              new_match.ref_start = i - old_layer[j] + 1;
              new_match.ref_end = i;
              new_match.str_start = j - old_layer[j] + 2;
              new_match.str_end = j + 1;
              matches.push_back (new_match);
            }
        }
      FREE_ARRAY (old_layer);
      FREE_ARRAY (new_layer);
      my_system->data[index].str.clear ();
      my_system->add_matches (index, matches);
    }
  pthread_barrier_wait (sync);

  for (size_t index = 0; index < size; ++index)
    {
      const string &str = my_system->data[index].str;
      if (str.empty ())
        continue;

      vector<matches_pair> matches;
      if (ref_tree)
        {
          find_matches (thread_index, thread_count, min_match_length, ref_tree,
                        ref, str, true, matches);
        }
      else
        {
          if (thread_index == 0)
            {
              str_tree = new prefix_tree (part_len);
              str_tree->parse (str);
            }
          pthread_barrier_wait (sync);

          str_tree->calculate_depth (thread_index, thread_count, sync);

          pthread_barrier_wait (sync);

          find_matches (thread_index, thread_count, min_match_length, str_tree,
                        str, ref, false, matches);
        }
      sort (matches.begin (), matches.end ());
      my_system->add_matches (index, matches);
      pthread_barrier_wait (sync);
      /// Calculation done
      if (thread_index + 1 == thread_count)
        /// TODO: Make parallel, merge sort needed to combine matches
        sort (my_system->data[index].matches.begin (), my_system->data[index].matches.end ());

      if (thread_index == 0)
        {
          FREE (str_tree);
        }
    }
  if (thread_index == 0)
    FREE (ref_tree);
}
