from handlers.BamHandler import BamHandler
from handlers.FastaHandler import FastaHandler
from collections import defaultdict
import sys


MAX_COVERAGE = 50
SEQUENCE_LENGTH_CUTOFF_FACTOR = 3   # exclude aligned segments that exceed window_size by this multiple
DEFAULT_MIN_MAP_QUALITY = 20

DELETE_CHAR = "*"
INSERT_CHAR = "_"


class PileupGenerator:
    def __init__(self, chromosome_name, start_position, end_position, ref_sequence, reads):
        # candidate finder includes end position, so should the reference sequence
        self.chromosome_name = chromosome_name
        self.start_position = start_position
        self.end_position = end_position
        self.window_size = end_position - start_position
        self.ref_sequence = ref_sequence
        self.max_coverage = MAX_COVERAGE

        self.reads = reads

        self.read_start_indices = dict()
        self.read_end_indices = dict()
        self.read_alignment_starts = dict()
        self.read_alignment_ends = dict()

        self.sequences = defaultdict(list)
        self.aligned_segments = defaultdict(list)
        self.qualities = defaultdict(list)
        self.cigars = defaultdict(list)

        self.read_insert_lengths = defaultdict(dict)
        self.positional_insert_lengths = dict()

        self.insert_offsets = defaultdict(int)

    def get_read_segments(self):
        for r,read in enumerate(self.reads):
            if read.mapping_quality >= DEFAULT_MIN_MAP_QUALITY and not read.is_secondary and not read.is_unmapped \
                    and not read.is_qcfail and not read.is_supplementary:

                if read.is_read2:
                    sys.stderr.write("WARNING: 'is_read2' flag found 'True' for read: %s" % read.query_name)
                    continue

                self.get_aligned_segment_from_read(read)

        return self.sequences

    def get_aligned_read_segments(self):
        for r,read in enumerate(self.reads):
            if read.mapping_quality >= DEFAULT_MIN_MAP_QUALITY and not read.is_secondary and not read.is_unmapped \
                    and not read.is_qcfail and not read.is_supplementary:

                read.query_name = read.query_name + '_1' if read.is_read1 else read.query_name + '_2'
                self.get_aligned_segment_from_read(read)

            if r == self.max_coverage:
                break

        self.add_insertion_padding()

        for read_id in self.aligned_segments:
            print(''.join(self.aligned_segments[read_id]))
            print(''.join(self.cigars[read_id]))

        return self.aligned_segments

    def add_insertion_padding(self):
        # print(self.start_position, self.end_position)
        # print(sorted(self.positional_insert_lengths.items()))

        for read_id in self.aligned_segments:
            for position in sorted(self.positional_insert_lengths):
                segment_index = position - self.start_position + self.insert_offsets[read_id]
                max_insert_length = self.positional_insert_lengths[position]

                if read_id in self.read_insert_lengths:
                    if position in sorted(self.read_insert_lengths[read_id]):
                        # read has insert here
                        read_insert_length = self.read_insert_lengths[read_id][position]
                        segment_index += read_insert_length
                        padding_length = max_insert_length - read_insert_length

                    else:
                        # read has insert... but not here
                        padding_length = max_insert_length
                else:
                    # read has no inserts
                    padding_length = max_insert_length

                # print(''.join(self.aligned_segments[read_id]))
                # print(''.join(self.cigars[read_id]))
                # print("position: ",position)
                # print("start: ", self.start_position)
                # print("padding: ", padding_length)
                # print("index: ", segment_index)

                sequence = self.aligned_segments[read_id]
                cigars = self.cigars[read_id]
                padded_sequence = sequence[:segment_index] + [INSERT_CHAR]*padding_length + sequence[segment_index:]
                padded_cigars = cigars[:segment_index] + [INSERT_CHAR]*padding_length + cigars[segment_index:]

                self.aligned_segments[read_id] = padded_sequence
                self.cigars[read_id] = padded_cigars

                self.insert_offsets[read_id] += max_insert_length

            # print(self.read_alignment_starts[read_id], self.read_alignment_ends[read_id])
            # print(self.start_position, self.end_position)
            # print(''.join(self.aligned_segments[read_id]))

    def update_positional_insert_lengths(self, position, length):
        if position not in self.positional_insert_lengths:
            self.positional_insert_lengths[position] = length
        else:
            previous_length = self.positional_insert_lengths[position]
            self.positional_insert_lengths[position] = max(previous_length, length)

    def get_aligned_segment_from_read(self, read):
        """
        :param read:
        :return:
        """

        read_alignment_start = read.reference_start
        # read_alignment_stop = self.get_read_stop_position(read)

        cigar_tuples = read.cigartuples
        read_sequence = read.query_sequence
        read_id = read.query_name
        # read_quality = read.query_qualities

        # read_index: index of read sequence
        # ref_index: index of reference sequence
        read_index = 0
        ref_index = 0
        found_valid_cigar = False
        completion_status = False

        if read_id in self.read_start_indices:
            print("WARNING: read_id hash conflict", read_id)

        for c,cigar in enumerate(cigar_tuples):
            cigar_code = cigar[0]
            length = cigar[1]

            # get the sequence segments that are affected by this operation
            # read_quality_segment = read_quality[read_index:read_index + length]
            read_sequence_segment = read_sequence[read_index:read_index + length]

            # skip parsing the first segment if it is not a match
            if cigar_code != 0 and found_valid_cigar is False:
                # only increment the read index if the non-match cigar code is INS or SOFTCLIP
                if cigar_code == 1 or cigar_code == 4:
                    read_index += length
                continue
            found_valid_cigar = True

            # send the cigar tuple to get data from this segment
            ref_index_increment, read_index_increment, completion_status = \
                self.parse_cigar_tuple(read_index=read_index,
                                       cigar_code=cigar_code,
                                       length=length,
                                       alignment_position=read_alignment_start + ref_index,
                                       read_segment=read_sequence_segment,
                                       read_id=read_id,
                                       completion_status=completion_status)

            # increase the read index iterator
            read_index += read_index_increment
            ref_index += ref_index_increment

            if completion_status or c == len(cigar_tuples)-1:
                start_index = self.read_start_indices[read_id]
                end_index = self.read_end_indices[read_id]

                segment_alignment_start = self.read_alignment_starts[read_id]
                segment_alignment_end = self.read_alignment_ends[read_id]

                # to simulate Paolo Carnevali's data, all reads should span the full region, match on start and end pos.
                if segment_alignment_start == self.start_position and segment_alignment_end == self.end_position:
                    sequence = read_sequence[start_index:end_index + 1]

                    if len(sequence) < SEQUENCE_LENGTH_CUTOFF_FACTOR*self.window_size:
                        self.sequences[read_id] = sequence

                else:
                    del self.aligned_segments[read_id]

                break

        return True

    def parse_delete(self, read_id, read_index, alignment_position, length):
        """
        Process a cigar operation that is a delete
        :param alignment_position: Alignment position
        :param length: Length of the delete
        :param ref_sequence: Reference sequence of delete
        :return:
        """
        # actual delete position starts one after the anchor
        read_complete = False
        index = read_index
        start = alignment_position
        stop = start + length

        for i in range(start, stop):
            in_left_bound = i >= self.start_position
            in_right_bound = i <= self.end_position

            # update start position
            if in_left_bound:
                if read_id not in self.read_start_indices:
                    self.read_start_indices[read_id] = index
                    self.read_alignment_starts[read_id] = i

            # update read
            if in_left_bound and in_right_bound:
                self.cigars[read_id].append("D")
                self.aligned_segments[read_id].append(DELETE_CHAR)

            # update end position
            if in_right_bound:
                self.read_end_indices[read_id] = index
                self.read_alignment_ends[read_id] = i
            else:
                read_complete = True
                break

            index += 1

        return read_complete

    def parse_insert(self, read_index, read_id, alignment_position, length, read_segment):
        """
        Process a cigar operation where there is an insert
        :param alignment_position: Position where the insert happened
        :param read_segment: The insert read sequence
        :return:
        """
        index = read_index
        read_complete = False
        start = alignment_position + 1
        stop = start + length

        segment_index = 0
        for i in range(start, stop):
            in_left_bound = i > self.start_position     # insert should not be first element in read
            in_right_bound = i <= self.end_position

            # update read
            if in_left_bound and in_right_bound:
                # Only update IFF this is not the first element in the alignment (should not start on insert)
                if read_id in self.read_start_indices:
                    character = read_segment[segment_index]
                    self.cigars[read_id].append("I")
                    self.aligned_segments[read_id].append(character)

                    if segment_index == 0:
                        position = start + segment_index - 1
                        self.read_insert_lengths[read_id][position] = length
                        self.update_positional_insert_lengths(position=position, length=length)

            # update end position
            if in_right_bound:
                self.read_end_indices[read_id] = index
                # self.read_alignment_ends[read_id] = i
            else:
                read_complete = True
                break

            index += 1
            segment_index += 1

        return read_complete

    def parse_match(self, read_index, read_id, alignment_position, length, read_segment):
        index = read_index
        read_complete = False
        start = alignment_position
        stop = start + length

        segment_index = 0
        for i in range(start, stop):
            in_left_bound = i >= self.start_position
            in_right_bound = i <= self.end_position

            # update start position
            if in_left_bound:
                if read_id not in self.read_start_indices:
                    self.read_start_indices[read_id] = index
                    self.read_alignment_starts[read_id] = i

            # update read
            if in_left_bound and in_right_bound:
                character = read_segment[segment_index]
                self.cigars[read_id].append("M")
                self.aligned_segments[read_id].append(character)

            # update end position
            if in_right_bound:
                self.read_end_indices[read_id] = index
                self.read_alignment_ends[read_id] = i
            else:
                read_complete = True
                break

            index += 1
            segment_index += 1

        return read_complete

    def parse_cigar_tuple(self, read_index, cigar_code, length, alignment_position, read_segment, read_id, completion_status):
        """
        Parse through a cigar operation to find possible candidate variant positions in the read
        :param cigar_code: Cigar operation code
        :param length: Length of the operation
        :param alignment_position: Alignment position corresponding to the reference
        :param read_segment: Read sequence from the bounds of this cigar operation
        :return:
        cigar key map based on operation.
        details: http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        0: "MATCH",
        1: "INSERT",
        2: "DELETE",
        3: "REFSKIP",
        4: "SOFTCLIP",
        5: "HARDCLIP",
        6: "PAD"
        """
        ref_index_increment = length
        read_index_increment = length

        # deal different kinds of operations
        if cigar_code == 0:
            # match
            completion_status = self.parse_match(read_index=read_index,
                                                 read_id=read_id,
                                                 alignment_position=alignment_position,
                                                 length=length,
                                                 read_segment=read_segment)

        elif cigar_code == 1:
            # insert
            completion_status = self.parse_insert(read_index=read_index,
                                                  read_id=read_id,
                                                  alignment_position=alignment_position,
                                                  length=length,
                                                  read_segment=read_segment)

            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            ref_index_increment = 0

        elif cigar_code == 2 or cigar_code == 3:
            # delete or ref_skip
            completion_status = self.parse_delete(read_index=read_index,
                                                  read_id=read_id,
                                                  alignment_position=alignment_position,
                                                  length=length)

            # alignment position is where the next alignment starts, for insert and delete this
            # position should be the anchor point hence we use a -1 to refer to the anchor point
            read_index_increment = 0

        elif cigar_code == 4:
            # soft clip
            ref_index_increment = 0
            # print("CIGAR CODE ERROR SC")

        elif cigar_code == 5:
            # hard clip
            ref_index_increment = 0
            read_index_increment = 0
            # print("CIGAR CODE ERROR HC")

        elif cigar_code == 6:
            # pad
            ref_index_increment = 0
            read_index_increment = 0
            # print("CIGAR CODE ERROR PAD")

        else:
            raise ("INVALID CIGAR CODE: %s" % cigar_code)

        return ref_index_increment, read_index_increment, completion_status