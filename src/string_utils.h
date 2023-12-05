/** @file string_utils.h
 *  @brief Header containing the StringUtils class
 *
 *  This contains the subroutines and eventually any 
 *  macros, constants, etc. needed for StringUtils class
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



#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include "constants.h"

namespace StringUtils
{
	
/** @brief It capitalizes a char_t type variable.
*
*  @param ch char_t type variable
*  @return Capitalized value
*/	
	Constants::char_t up_char(Constants::char_t ch);
	
	
/** @brief It capitalizes every char of a string.
*
*  @param src String
*  @return Capitalized string
*/	
	std::string toupper(const std::string &src);
	
	
/** @brief It converts a char_t type variable into lower case.
*
*  @param ch char_t type variable
*  @return Lower case value
*/	
	Constants::char_t down_char(Constants::char_t ch);
	
	
/** @brief It converts every char of a string into lower case.
*
*  @param src String
*  @return Lower case string
*/		
	std::string tolower(const std::string &src);
	
	
/** @brief It determines a string is a numeric number or not.
*
*  @param src String
*  @return True or False
*/	
	bool is_numeric(const std::string& str);


/** @brief It splits a string by a char delimeter.
*
*  @param s String
*  @param delim Char delimeter
*  @param elems String vector
*  @return Splited string vector
*/
	Constants::string_vector &split(const std::string &s, char delim, Constants::string_vector &elems);
	
	
/** @brief It splits a string by a char delimeter.
*
*  @param s String
*  @param delim Char delimeter
*  @return Splited string vector
*/	
	Constants::string_vector split(const std::string &s, char delim);
	
	
/** @brief It converts every element of a string vector into an integer vector.
*
*  @param vs String vector
*  @return Integer vector
*/	
	std::vector<int> vecstr_to_vecint(std::vector<std::string> vs);
	
	
/** @brief It converts every element of a string vector into an floating point vector.
*
*  @param vs String vector
*  @return Floating point vector
*/		
	template <typename T>
	std::vector<T> vecstr_to_vecflt(Constants::string_vector vs);


/** @brief It converts a integer to string.
*
*  @param i Integer number
*  @return String
*/	
	std::string itoa(int i);
	
	
/** @brief It converts a integer to string.
*
*  @param num Integer number
*  @return String
*/	
	std::string itos(int num);


/** @brief It trims left side of a string.
*
*  @param s Input string
*  @return Left trimmed string
*/
	static inline std::string& ltrim(std::string &s)
	{
		s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		return s;
	}


/** @brief It trims right side of a string.
*
*  @param s Input string
*  @return Right trimmed string
*/
	static inline std::string& rtrim(std::string &s)
	{
		s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
		return s;
	}


/** @brief It trims a string from both side.
*
*  @param s Input string
*  @return Trimmed string
*/
	static inline std::string& trim(std::string &s)
	{
		return ltrim(rtrim(s));
	}


	Constants::char_t up_char(Constants::char_t ch)
	{
		return std::use_facet<std::ctype<Constants::char_t>>(std::locale()).toupper(ch);
	}


	std::string toupper(const std::string &src)
	{
		std::string result;
		std::transform(src.begin(), src.end(), std::back_inserter(result), up_char);
		return result;
	}


	Constants::char_t down_char(Constants::char_t ch)
	{
		return std::use_facet<std::ctype<Constants::char_t>>(std::locale()).tolower(ch);
	}


	std::string tolower(const std::string &src)
	{
		std::string result;
		std::transform(src.begin(), src.end(), std::back_inserter(result), down_char);
		return result;
	}


	bool is_numeric(const std::string& str)
	{
		if (str.size() == 0)
		return false;

		std::stringstream conv;
		double tmp;
		conv << str;
		conv >> tmp;
		return conv.eof();
	}


	Constants::string_vector& split(const std::string &s, char delim, Constants::string_vector &elems)
	{
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim))
		{
			elems.push_back(item);
		}
		return elems;
	}


	Constants::string_vector split(const std::string &s, char delim)
	{
		std::vector<std::string> elems;
		return split(s, delim, elems);
	}


	std::vector<int> vecstr_to_vecint(Constants::string_vector vs)
	{
		std::vector<int> ret;
		for (Constants::string_vector::iterator it = vs.begin(); it != vs.end(); ++it)
		{
			std::istringstream iss(*it);
			int temp;
			iss >> temp;
			ret.push_back(temp);
		}
		return ret;
	}


	template <typename T>
	std::vector<T> vecstr_to_vecflt(Constants::string_vector vs)
	{
		std::vector<T> ret;
		for (Constants::string_vector::iterator it = vs.begin(); it != vs.end(); ++it)
		{
			std::istringstream iss(*it);
			T temp;
			iss >> temp;
			ret.push_back(temp);
		}
		return ret;
	}


	std::string itoa(int i)
	{
		return std::to_string(i);
	}


	std::string itos(int num)
	{
		std::ostringstream ss;
		ss << num;
		return ss.str();
	}
}

#endif
