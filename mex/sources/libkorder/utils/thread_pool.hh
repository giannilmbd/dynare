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

#ifndef THREAD_POOL_HH
#define THREAD_POOL_HH

#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace thread_pool
{
class job
{
public:
  virtual ~job() = default;
  // Run the job. The mutex is shared between all jobs of the same group.
  virtual void operator()(std::mutex& mut) = 0;
};

// Spawn the threads in the pool. May be called several times (in which case the method ensures that
// the right number of threads is present).
void initialize(); // Use half the number of logical cores as number of threads
void initialize(int thread_number);

[[nodiscard]] int get_thread_number();

using job_group_t = std::vector<std::unique_ptr<job>>;

// Run a group of jobs, and only returns when they are all completed
void run(const job_group_t& job_group);
}

#endif
