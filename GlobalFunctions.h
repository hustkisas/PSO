#ifndef GLOBALFUNCTIONS_H_HEADER_INCLUDED_BECEF599
#define GLOBALFUNCTIONS_H_HEADER_INCLUDED_BECEF599
#include "Global.h"

static std::vector<std::string> ReadSectionsFromFile(std::string& fileName, std::string& front_delim, std::string& end_delim)
{
	std::ifstream   inputFile( fileName.c_str()); 
	std::string   fileData((std::istreambuf_iterator <char> (inputFile)), std::istreambuf_iterator <char> ()); 
	std::vector<std::string> sections;
	std::string delim = front_delim+end_delim;
	size_t position = 0;
	size_t pos_section_begin = -1;
	size_t pos_section_end = -1;
	int cnt_front_delim = 0;
	// 	std::string subString = fileData;
	while (position != std::string::npos || position == fileData.length())
	{
		position = fileData.find_first_of(delim, position+1);
		if (position == std::string::npos)
		{
			break;
		}
		std::string chr = fileData.substr(position, 1);
		if (chr == front_delim)
		{
			if (cnt_front_delim ==0)
				pos_section_begin = position;

			cnt_front_delim++;
		}
		else if (chr == end_delim)
		{
			cnt_front_delim--;
			if (cnt_front_delim == 0)
				pos_section_end = position;
		}


		if (pos_section_begin != -1 && pos_section_end != -1)
		{
			std::string section = fileData.substr(pos_section_begin, pos_section_end - pos_section_begin);
			section = section.substr(1, section.length() - 1);
			sections.push_back(section);
			position = pos_section_end;
			pos_section_begin = -1;
			pos_section_end = -1;
		}
	}
	return sections;
}

static std::vector<std::string> ReadSectionsFromString(std::string& strContent, std::string& front_delim, std::string& end_delim)
{
	std::vector<std::string> sections;
	std::string delim = front_delim+end_delim;
	size_t position = strContent.find_first_of(delim);
	size_t pos_section_begin = -1;
	size_t pos_section_end = -1;
	int cnt_front_delim = 0;
	// 	std::string subString = fileData;
	while (position != std::string::npos || position == strContent.length())
	{
		if (position == std::string::npos)
		{
			break;
		}
		std::string chr = strContent.substr(position, 1);
		if (chr == front_delim)
		{
			if (cnt_front_delim ==0)
				pos_section_begin = position;

			cnt_front_delim++;
		}
		else if (chr == end_delim)
		{
			cnt_front_delim--;
			if (cnt_front_delim == 0)
				pos_section_end = position;
		}


		if (pos_section_begin != -1 && pos_section_end != -1)
		{
			std::string section = strContent.substr(pos_section_begin, pos_section_end - pos_section_begin);
			section = section.substr(1, section.length() - 1);
			sections.push_back(section);
			position = pos_section_end;
			pos_section_begin = -1;
			pos_section_end = -1;
		}

		position = strContent.find_first_of(delim, position+1);

	}
	return sections;
}



static void   Search(int sizeLimit, int   m,   int   n,   int   depth,   VEC_INT&   mark,   VEC_INT&   L,   std::list<VEC_INT>&   result)  
{  
	if (result.size() >= sizeLimit)
		return;

	if(   depth   ==   n   )   {  
		result.push_back(L);  
	}   
	else   {  
		int   begin;  
		if(   depth   ==   0   )   {  
			begin   =   0;  
		}   
		else   {  
			begin   =   L[depth-1];  
		}  
		for(   int   i   =   begin;   i   <   m;   i++   )    
			if(   mark[i]   ==   0   )   {  
				mark[i]   =   1;  
				L[depth]   =   i  ; // +   1;  
				Search(sizeLimit,   m,   n,   depth+1,   mark,   L,   result   );  
				mark[i]   =   0;  
			}  
	}  
}  

static void   Combination(int sizeLimit, int   m,   int   n,   std::list<VEC_INT>&   result)  
{                          
	if(   m   <   n   )   return;  

	VEC_INT   mark(m);                       //   用来标记某个数字是否用过了  

	for(   int   i=0;   i   <   m;   i++   )   {   //   初始化mark  
		mark.push_back(0);  
	}  

	VEC_INT   comb(n);               //   存储一组组合数      
	comb.resize(n);  

	Search(sizeLimit, m,   n,   0,   mark,   comb,   result);  
}  

#endif

