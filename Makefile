
target := bhc npg

all: $(addsuffix .html,$(target)) 

$(addsuffix .html,$(target)) : %.html : %.jemdoc imgRecSrc.conf MENU
	jemdoc -c imgRecSrc.conf $<

