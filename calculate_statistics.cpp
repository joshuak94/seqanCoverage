#include <seqan3/io/sam_file/input.hpp>     // SAM/BAM support
#include <seqan3/argument_parser/all.hpp>   // For the argument parsing
#include <seqan3/core/debug_stream.hpp>

#include <seqan/bam_io.h>

#include <filesystem>
#include <iostream>
#include <sstream>
#include <chrono>
#include <ctime>
#include <cstring>

using namespace seqan3::literals;

#define RUN(x, y) {startTimeMessage(y);x;endTimeMessage(y);}

time_t _t1, _t2;
std::chrono::time_point<std::chrono::system_clock> _m1, _m2;

void printTimeMessage(std::string msg)
{
    time_t t = time(0);
    struct tm* now = localtime(&t);
    seqan3::debug_stream << "[ " << std::put_time(now, "%d-%m-%Y %H:%M:%S") << "] " << msg << '\n';
}

void startTimeMessage(std::string msg)
{
    time(&_t1);
    _m1 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
    printTimeMessage("[START] " + msg);
}

void endTimeMessage(std::string msg)
{
    using std::chrono::operator""ms;
    using std::chrono::operator""us;
    time(&_t2);
    _m2 = std::chrono::time_point_cast<std::chrono::microseconds>(std::chrono::system_clock::now());
    auto difference = std::chrono::duration_cast<std::chrono::microseconds>(_m2 - _m1);
    if (difference < 1000us)
    {
        printTimeMessage("[END] " + msg + " (" + std::to_string(difference.count()) + " microseconds.)");
    }
    else if (difference >= 1000us && difference < 1000ms)
    {
        printTimeMessage("[END] " + msg + " (" +
                         std::to_string((std::chrono::duration_cast<std::chrono::milliseconds>(difference)).count()) +
                         " milliseconds.)");
    }
    else
    {
        printTimeMessage("[END] " + msg + " (" + std::to_string((int)difftime(_t2,_t1)) + " seconds.)");
    }
}

struct CmdOptions
{
    std::filesystem::path input_path{};
    uint32_t num_threads{1};
};

void initialize_parser(seqan3::argument_parser & parser, CmdOptions & options)
{
    parser.add_option(options.input_path, 'i', "input",
                      "The path to the input BAM/SAM file. Must be sorted!",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"sam", "bam"}});
    parser.add_option(options.num_threads, 't', "threads",
		      "The number of threads to use for decompression.");
}

std::vector<std::uint64_t> parse_sam_file_default(std::filesystem::path const & input_path)
{
    seqan3::sam_file_input input_file{input_path};
    std::vector<std::uint64_t> ref_offset_vector{};
    for (auto & rec : input_file)
    {
        ref_offset_vector.push_back(rec.reference_position().value_or(0));
    }

    return ref_offset_vector;
}

std::vector<std::uint64_t> parse_sam_file_minimal(std::filesystem::path const & input_path)
{
    seqan3::sam_file_input input_file{input_path, seqan3::fields<seqan3::field::ref_offset>{}};
    std::vector<uint64_t> ref_offset_vector{};
    for (auto & rec : input_file)
    {
        ref_offset_vector.push_back(rec.reference_position().value_or(0));
    }

    return ref_offset_vector;
}

std::vector<std::uint64_t> parse_sam_file_seqan2(std::filesystem::path const & input_path)
{
    seqan::BamFileIn input_file;
    std::vector<uint64_t> ref_offset_vector{};
    if (!open(input_file, seqan::toCString(input_path.string())))
    {
        std::cerr << "ERROR: Could not open " << input_path << std::endl;
        return ref_offset_vector;
    }

    try
    {
        // Copy header.
        seqan::BamHeader header;
        seqan::readHeader(header, input_file);

        // Copy records.
        seqan::BamAlignmentRecord record;
        while (!seqan::atEnd(input_file))
        {
            seqan::readRecord(record, input_file);
            ref_offset_vector.push_back(record.beginPos);
        }
    }
    catch (seqan::Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
    }
    return ref_offset_vector;
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

    seqan3::contrib::bgzf_thread_count = options.num_threads;
    RUN(parse_sam_file_seqan2(options.input_path), "Seqan2 parsing.");
    RUN(parse_sam_file_minimal(options.input_path), "Minimal parsing.");
    RUN(parse_sam_file_default(options.input_path), "Default parsing.");

}
