/** @file supertimer.h
 *  @brief Header containing the SuperTimer class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for SuperTimer class
 *
 *  @author Mario Morales Hernandez
 *  @author Md Bulbul Sharif
 *  @author Tigstu T. Dullo
 *  @author Sudershan Gangrade
 *  @author Alfred Kalyanapu
 *  @author Sheikh Ghafoor
 *  @author Shih-Chieh Kao
 *  @author Katherine J. Evans
 *  @bug No known bugs.
 */



#ifndef SUPERTIMER_H
#define SUPERTIMER_H

#include "constants.h"

namespace SuperTimer
{
	struct ci_less	/**< Structure to compare to string. */
	{
		struct nocase_compare	/**< Structure to compare to char. */
		{
			
/** @brief It compares to char. If first char is less than second char, it returns true.
*
*  @param c1 First char
*  @param c2 Second char
*  @return True or False
*/
			bool operator() (const unsigned char& c1, const unsigned char& c2) const
			{
				return tolower(c1) < tolower(c2);
			}
		};
		
		
/** @brief It compares to string. If first string is less than second string, it returns true. Compare happens one by one char.
*
*  @param s1 First string
*  @param s2 Second string
*  @return True or False
*/		
		bool operator() (const std::string & s1, const std::string & s2) const
		{
			return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end(), nocase_compare());
		}
	};


	class super_timer	/**< Custom timer class to compute time for different operation. */
	{
	public:
		
/** @brief Constructor
*
*/	
		super_timer();


/** @brief It starts a custom timer.
*
*  @param category Timer name
*/
		void start(std::string category);
		
		
/** @brief It stops a custom timer.
*
*  @param category Timer name
*/		
		void stop(std::string category);

/** @brief It restarts a custom timer (set to zero)
*
*  @param category Timer name
*/
		void restart(std::string category);	
		
/** @brief It resets all timers.
*
*/		
		void reset();


/** @brief It calculates total time of every timer.
*
*  @return Time value
*/
		double get_total_time();
		
		
/** @brief It calculates time value of a specific timer.
*
*  @param category Timer name
*  @return Time value
*/		
		double get_custom_time(std::string category);


/** @brief It calculates current date.
*
*  @return Current date
*/
		std::string get_current_date();
		
		
/** @brief It calculates host name.
*
*  @return Host name
*/		
		std::string get_hostname();


/** @brief It helps to add a new timer.
*
*  @param category Timer name
*  @return Timer id
*/
		int add_new_timer(std::string category);


	private:
		int units_;	/**< Time units of a timer */
		double unit_factor_;	/**< Time unit factor of a timer */
		std::vector<timeval> all_timers_;	/**< Vector that contains all timer information */
		std::map<std::string, int, ci_less> timecats_;	/**< Timer name, size mapping */
		std::vector<Constants::ull> times_;	/**< Vector to contains timers time value */


/** @brief It sets time of a timer.
*
*  @param category Timer name
*  @param time Timer time value
*/
		void upd_time_(std::string category, Constants::ull time);
		
		
/** @brief It calculates index id of a timer.
*
*  @param category Timer name
*  @return Timer id
*/		
		int get_cat_index_(std::string category);
		
		
/** @brief It initializes timer object.
*
*/		
		void init_();
		
		
/** @brief It updates time of a timer.
*
*  @param timer Timer id
*  @param time Timer time value
*/		
		void upd_time_(int timer, Constants::ull time);
		
		
/** @brief It converts time value to floating point number.
*
*  @param time Timer time value
*  @return Time value
*/		
		double convert_(Constants::ull time);
	};


	super_timer::super_timer() : units_(TIMER_SECS)
	{
		init_();
	}


	void super_timer::init_()
	{
		unit_factor_ = 1;
		if(units_ == TIMER_SECS)
		{
			unit_factor_ = 0.000001;
		}
		else if (units_ == TIMER_NSECS)
		{
			unit_factor_ = 1000;
		}
	}


	int super_timer::add_new_timer(std::string category)
	{
		int ncats = timecats_.size();
		std::map<std::string, int, ci_less>::iterator it = timecats_.find(category);
		if (it == timecats_.end())
		{
			timeval ntval;
			Constants::ull newcat_time = 0;

			timecats_.insert(std::pair<std::string, int>(category, ncats));
			times_.push_back(newcat_time);
			all_timers_.push_back(ntval);
		}

		return timecats_.size() - 1;
	}


	int super_timer::get_cat_index_(std::string category)
	{
		int catidx;
		std::map<std::string, int, ci_less>::iterator it = timecats_.find(category);

		if (it != timecats_.end())
		{
			catidx = it->second;
		}
		else
		{
			add_new_timer(category);
			catidx = add_new_timer(category);
		}

		return catidx;
	}


	void super_timer::upd_time_(int catid, Constants::ull time)
	{
		if (catid < (int)times_.size())
		{
			times_[catid] += time;
		}
	}


	void super_timer::upd_time_(std::string category, Constants::ull time)
	{
		times_[get_cat_index_(category)] += time;
	}


	double super_timer::convert_(Constants::ull time)
	{
		return (double)time * unit_factor_;
	}


	void super_timer::reset()
	{
		all_timers_.clear();
		times_.clear();

		init_();
	}


	void super_timer::start(std::string category)
	{
		timeval* t = &all_timers_[get_cat_index_(category)];

		gettimeofday(t, NULL);
	}


	void super_timer::stop(std::string category)
	{
		timeval t1, t2;
		int catidx = get_cat_index_(category);
		gettimeofday(&t2, NULL);
		t1 = all_timers_[catidx];

		Constants::ull time_usecs = ((((Constants::ull) t2.tv_sec * 1000000) + t2.tv_usec) - (((Constants::ull) t1.tv_sec * 1000000) + t1.tv_usec));

		upd_time_(category, time_usecs);
	}

	void super_timer::restart(std::string category)
	{
		times_[get_cat_index_(category)] = 0.0;
	}



	double super_timer::get_custom_time(std::string category)
	{
		return convert_(times_[get_cat_index_(category)]);
	}


	double super_timer::get_total_time()
	{
		double total = 0.0;
		std::map<std::string, int>::iterator it = timecats_.begin();

		for (it = timecats_.begin(); it != timecats_.end(); it++)
		{
			total += get_custom_time(it->first);
		}

		return total;
	}

	std::string super_timer::get_current_date()
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[24];

		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(buffer, 24, "%D", timeinfo);

		return std::string(buffer);
	}


	std::string super_timer::get_hostname()
	{
		char hostname[1024];
		hostname[1023] = '\0';
		gethostname(hostname, 1023);

		return std::string(hostname);
	}
}

#endif
