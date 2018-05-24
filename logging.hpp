#pragma once

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;

typedef sinks::synchronous_sink< sinks::text_file_backend > file_sink;

namespace ds2i {


#define DS2I_LOG BOOST_LOG_TRIVIAL(info)

void init_logging()
{
    boost::log::core::get()->remove_all_sinks();
    boost::log::core::get()->add_global_attribute("TimeStamp", boost::log::attributes::local_clock());
    logging::add_console_log(std::cout,keywords::format = "[%TimeStamp%]: %Message%");
}

void start_logging_to_file(std::string filename)
{
    boost::log::core::get()->remove_all_sinks();
    boost::log::core::get()->add_global_attribute("TimeStamp", boost::log::attributes::local_clock());
    logging::add_console_log(std::cout,keywords::format = "[%TimeStamp%]: %Message%");
    logging::add_file_log
	    (
	        keywords::file_name = filename,
	        keywords::format = "[%TimeStamp%]: %Message%",
		keywords::auto_flush = true
	    );
}

void stop_logging_to_file()
{
    boost::log::core::get()->remove_all_sinks();
    boost::log::core::get()->add_global_attribute("TimeStamp", boost::log::attributes::local_clock());
    logging::add_console_log(std::cout,keywords::format = "[%TimeStamp%]: %Message%");
}

}
