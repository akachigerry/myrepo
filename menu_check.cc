#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

void menuList (const std::string *array, std::size_t len)  	//when passed to a function, an array 'decays' to a pointer, 
																	//hence cannot use its size member. To use array size in function, 
																		//its length must be supplied as parameter in function definition
																			//and as argument when function is used.

{   
	for (size_t i=0; i < len; i++)			//this step moves pointer from first item in array to last item. 
										//size_t is type for representing byte sizes
		std::cout << array[i] << '\n';
}
	
int main() {
	std::vector<std::string> menu = {"Bread_Bacon", "Bread_Omlette", "Chicken_Biryani", "Chicken_teryaki", "Lamb_Biryani", "Chicken_korma"};
	std::cout << "Heres your menu for today" << std::endl;
	menuList(menu.data(), menu.size());			//menu.data() is reference to first item of menu array.
	std::cout << "Please choose from above menu" << "\n";
	std::string menuChoice;
	std::cout << ">>";
	std::cin >> menuChoice;
	
	/*
	Check if menuChoice is in list of menu
	*/
	while (1) {
		std::vector<std::string>::const_iterator it = std::find(menu.begin(), menu.end(), menuChoice);
		if (it != menu.end()) {	
			std::cout << "You selected" << menuChoice << std::endl;
			break;
		} else {
			std::cout << "You have not selected from our menu list. Please select from the the following menu:" << "\n";
			menuList(menu.data(), menu.size());
			std::cout << ">>";
			std::cin >> menuChoice;
		}
	}
	return 0;
}
