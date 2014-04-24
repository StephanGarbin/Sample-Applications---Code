#pragma once

#ifndef DROGON_DROGON_STATUS_H
#define DROGON_DROGON_STATUS_H

#include <tbb\tick_count.h>
#include <iostream>

using namespace tbb;

namespace drogon {


	template<class T>
	class drogon_status
	{
	public:
		drogon_status(void)
		{
		}
		drogon_status(bool was_success, tick_count start, tick_count end)
		{
			this->was_success = was_success;
			this->start = start;
			this->end = end;
		}
		//for iterative solvers
		drogon_status(bool was_success, tick_count start, tick_count end, size_t iterations, T inf_norm_residual, T tolerance)
		{
			this->was_success = was_success;
			this->start = start;
			this->end = end;
			this->iterations = iterations;
			this->inf_norm_residual = inf_norm_residual;
			this->tolerance = tolerance;
		}

		~drogon_status()
		{
		}

		void print_cg_results()
		{
			if (was_success)
			{
				cout << "\n" << "Solution converged. [ITERATIONS | TIME | INF-NORM RESIDUAL | TOL]" << "\n";
				cout << "[" << iterations << " | " << (end - start).seconds() << "s" << " | " << inf_norm_residual << " | " << tolerance << "]";
				cout << "\n";
			}
			else
			{
				cout << "\n" << "Solution did not converge. [ITERATIONS | TIME | INF-NORM RESIDUAL | TOL]" << "\n";
				cout << "[" << iterations << " | " << (end - start).seconds() << "s" << " | " << inf_norm_residual << " | " << tolerance << "]" << "\n";
				cout << "\n";
			}
		}

		bool get_success()
		{
			return was_success;
		}

		void set_succes(bool success)
		{
			this->was_success = success;
		}

		tick_count get_start()
		{
			return start;
		}

		tick_count get_end()
		{
			return end;
		}

		size_t get_iterations()
		{
			return iterations;
		}

		T get_infnorm_residual()
		{
			return inf_norm_residual;
		}

		T get_tolerance()
		{
			return tolerance;
		}

		void set_start(tick_count start)
		{
			this->start = start;
		}

		void set_end(tick_count end)
		{
			this->end = end;
		}

		void set_iterations(size_t iterations)
		{
			this->iterations = iterations;
		}

		void set_infnorm_residual(T inf_norm)
		{
			this->inf_norm_residual = inf_norm;
		}

		void set_tolerance(T tolerance)
		{
			this->tolerance = tolerance;
		}


	private:
		bool was_success;
		tick_count start;
		tick_count end;

		//for iterative solvers
		size_t iterations;
		T inf_norm_residual;
		T tolerance;
	};



} /*namespace drogon*/



#endif /*DROGON_DROGON_STATUS_H*/