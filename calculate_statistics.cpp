#include <seqan3/io/sam_file/input.hpp>     // SAM/BAM support
#include <seqan3/argument_parser/all.hpp>   // For the argument parsing
#include <seqan3/core/debug_stream.hpp>

#include <filesystem>
#include <iostream>
#include <sstream>

using seqan3::operator""_cigar_op;

struct CmdOptions
{
    std::filesystem::path input_path{};
    std::filesystem::path output_path{};
};

template<typename stream_type>
//!cond
    requires seqan3::output_stream<stream_type>
//!endcond
void write_output(std::string chrName, std::vector<uint32_t> & chrCov, stream_type & output_stream)
{
    int pos{1};
    std::stringstream buffer{};
    for (auto i : chrCov)
    {
        if (i > 0)
        {
            buffer << chrName << '\t' << pos << '\t' << i << '\n';
        }
        ++pos;
    }
    output_stream << buffer.rdbuf();
}

void initialize_parser(seqan3::argument_parser & parser, CmdOptions & options)
{
    parser.add_option(options.input_path, 'i', "input",
                      "The path to the input BAM/SAM file. Must be sorted!",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"sam", "bam"}});
    parser.add_option(options.output_path, 'o', "output",
                      "The path and filename for the desired output file. Default is to output to standard out.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {}});
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"SeqAnCoverage", argc, argv};
    CmdOptions options{};
    initialize_parser(parser, options);
    try
    {
         parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "[ERROR] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    // Open input alignment file
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::cigar,
                                     seqan3::field::flag,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::mapq>;

    seqan3::sam_file_input alignment_file{options.input_path, my_fields{}};
    std::ofstream output_stream{};

    if (!options.output_path.empty())
    {
        output_stream.open(options.output_path.c_str());
    }

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
                if (options.output_path.empty())
                {
                    write_output(alignment_file.header().ref_ids()[curChr], chrCov, std::cout);
                }
                else
                {
                    write_output(alignment_file.header().ref_ids()[curChr], chrCov, output_stream);
                }
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
        for (size_t i = 0; i < match_len.size(); ++i)
        {
            std::transform(std::begin(chrCov) + pos + match_pos[i], std::begin(chrCov) + pos + match_pos[i] + match_len[i], std::begin(chrCov) + pos + match_pos[i],[](auto x){return x+1;});
        }
    }

    // Write last chromosome to file.
    if (options.output_path.empty())
    {
        write_output(alignment_file.header().ref_ids()[curChr], chrCov, std::cout);
    }
    else
    {
        write_output(alignment_file.header().ref_ids()[curChr], chrCov, output_stream);
        output_stream.close();
    }
}
