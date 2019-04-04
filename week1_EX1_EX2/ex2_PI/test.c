#include <stdio.h>

	int five() {
		return 5;
	}
int main() {

	for (int i = 0; i < 10; i++) {
		printf("%d ", five());
	}

	return 0;
}