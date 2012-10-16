
test:
	-cd Tests && ./Test.sh

clean:
	-@(find . -name "*.pyc" | xargs rm)

