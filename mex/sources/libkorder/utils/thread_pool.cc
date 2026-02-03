/*
 * Copyright © 2025-2026 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cassert>
#include <condition_variable>
#include <map>

#include "thread_pool.hh"

namespace thread_pool
{
namespace // Use an unnamed namespace to ensure that its contents cannot be accessed from
          // elsewhere
{
enum class status
{
  waiting,
  ongoing,
  done
};

/* NB: The threads must be listed *after* all the other variables, so that when MATLAB/Octave
 exits, the threads are destroyed *before* those variables (since the former rely on the latter
 for their functioning). */
std::mutex global_mut;                   // Protects the modification of the variables below
std::condition_variable_any thread_cv;   // Threads wait on this one
std::condition_variable_any director_cv; // The director waits on this one
std::map<job*, status> job_group_status;
std::mutex job_group_mut; // Shared by jobs in the job group

std::vector<std::jthread> threads;
} // End of unnamed namespace

void
initialize()
{
  initialize(std::max(1, static_cast<int>(std::thread::hardware_concurrency()) / 2));
}

void
initialize(int thread_number)
{
  assert(thread_number > 0);

  // Handle the case where initialize() is called several times
  // Happens e.g. in successive MEX calls
  if (static_cast<int>(threads.size()) == thread_number)
    return;
  else
    threads.clear();

  for (int i {0}; i < thread_number; i++)
    /* Passing the stop_token by const reference is ok (and makes clang-tidy happier),
       since the std::jthread constructor calls the lambda with the return argument of the
       get_stop_token() method, which returns a stop_token by value; hence there is no lifetime
       issue. See:
       https://stackoverflow.com/questions/72990607/const-stdstop-token-or-just-stdstop-token-as-parameter-for-thread-funct
     */
    threads.emplace_back([](const std::stop_token& stoken) {
      std::unique_lock lk {global_mut};
      job* selected_job;
      status* selected_job_status;

      auto pick_job = [&selected_job, &selected_job_status] {
        for (auto& [job, job_status] : job_group_status)
          if (job_status == status::waiting)
            {
              selected_job = job;
              selected_job_status = &job_status;
              job_status = status::ongoing;
              return true;
            }
        return false;
      };

      while (!stoken.stop_requested())
        if (thread_cv.wait(lk, stoken, pick_job))
          {
            lk.unlock();
            selected_job->operator()(job_group_mut);
            lk.lock();

            *selected_job_status = status::done;

            director_cv.notify_one();
          }
    });
}

int
get_thread_number()
{
  return threads.size();
}

void
run(const job_group_t& job_group)
{
  // Ensure that initialize() has been called
  assert(!threads.empty());
  // Ensure that run() is only called from the main thread, not from a parallel job
  assert(job_group_status.empty());

  std::unique_lock lk {global_mut};

  for (auto& job : job_group)
    job_group_status.try_emplace(job.get(), status::waiting);

  thread_cv.notify_all();

  director_cv.wait(lk, []() {
    return std::ranges::all_of(job_group_status, [](auto it) { return it.second == status::done; });
  });

  job_group_status.clear();
}
}
