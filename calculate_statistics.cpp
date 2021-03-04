#include <seqan3/io/sam_file/input.hpp>   // SAM/BAM support
#include <seqan3/core/debug_stream.hpp>

#include <filesystem>

using seqan3::operator""_cigar_op;

template<typename stream_type>
//!cond
    requires seqan3::output_stream<stream_type>
//!endcond
void write_output(std::string chrName, std::vector<uint32_t> & chrCov, stream_type & out_stream)
{
    int pos{1};
    for (auto i : chrCov)
    {
        if (i > 0)
        {
            out_stream << chrName << "\t" << pos << "\t" << i << std::endl;
        }
        pos++;
    }
}

int main([[ maybe_unused ]]int argc, [[ maybe_unused ]]char const ** argv)
{
    std::filesystem::path alignment_path = "/project/zvi/Joshua_Kim/VistaEnhancerValidation/PoolSequencing/P-PB.bam";
    std::filesystem::path output_file_path = "/home/kim_j/development/b_scratch/output.txt";
    std::ofstream out_file{output_file_path.c_str()};

    // Open input alignment file
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::cigar,
                                     seqan3::field::flag,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::mapq>;

    seqan3::sam_file_input alignment_file{alignment_path, my_fields{}};

    int32_t curChr{-1};
    std::vector<uint32_t> chrCov{};

    for (auto & rec : alignment_file)
    {
        std::string query_name                  = seqan3::get<seqan3::field::id>(rec);
        auto cigar_vector                       = seqan3::get<seqan3::field::cigar>(rec);                   // 6: CIGAR
        const seqan3::sam_flag flag             = seqan3::get<seqan3::field::flag>(rec);                    // 2: FLAG
        const auto ref_id                    = seqan3::get<seqan3::field::ref_id>(rec).value_or(-1);     // 3: RNAME
        const auto pos                       = seqan3::get<seqan3::field::ref_offset>(rec).value_or(-1);  // 4: POS
        std::vector<int32_t> match_len{};
        std::vector<int32_t> match_pos{};
        uint32_t skip{0};

        // Skip incorrect reads.
        if ((static_cast<uint16_t>(flag) & 0x0004) == 0x0004 || ref_id == -1 || pos == -1) continue;
        // Write to output before starting the next chromosome.
        if (curChr != ref_id)
        {
            if (curChr != -1)
            {
                write_output(alignment_file.header().ref_ids()[curChr], chrCov, out_file);
            }
            curChr = ref_id;
            chrCov.assign(std::get<0>(alignment_file.header().ref_id_info[curChr]), 0);
        }

        // Calculate covered positions from CIGAR string.
        for (seqan3::cigar & pair : cigar_vector)
        {
            using seqan3::get;
            int32_t length = get<0>(pair);
            seqan3::cigar_op operation = get<1>(pair);
            if (operation == 'S'_cigar_op || operation == 'H'_cigar_op || operation == 'I'_cigar_op || operation == 'P'_cigar_op) continue;
            if (operation == 'M'_cigar_op || operation == '='_cigar_op || operation == 'X'_cigar_op)
            {
                match_len.push_back(length);
                match_pos.push_back(skip);
            }
            // Okay for D and N
            skip += length;
        }

        // Add coverage to vector.
        for (size_t i = 0; i < match_len.size(); i++)
        {
            std::transform(std::begin(chrCov) + pos + match_pos[i], std::begin(chrCov) + pos + match_pos[i] + match_len[i], std::begin(chrCov) + pos + match_pos[i],[](auto x){return x+1;});
        }
    }
    write_output(alignment_file.header().ref_ids()[curChr], chrCov, out_file);
    out_file.close();

}
