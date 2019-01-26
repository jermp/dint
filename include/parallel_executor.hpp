#define BOOST_THREAD_VERSION 4
#define BOOST_THREAD_PROVIDES_EXECUTORS

#include <boost/config.hpp>
#include <boost/thread/executors/basic_thread_pool.hpp>
#include <boost/thread/experimental/parallel/v2/task_region.hpp>

typedef boost::executors::basic_thread_pool executor_type;
typedef boost::experimental::parallel::v2::task_region_handle_gen<executor_type>
    task_region_handle;
using boost::experimental::parallel::v2::task_region;

namespace ds2i {

    struct parallel_executor {
        parallel_executor(size_t num_threads = std::thread::hardware_concurrency()) {
            executor.reset(new executor_type(num_threads));
        }

        std::unique_ptr<executor_type> executor;
    };

}
