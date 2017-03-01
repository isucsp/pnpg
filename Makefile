
target := pnpg

all: $(addsuffix .html,$(target)) 

FLAGS := -c pnpg.conf

$(addsuffix .html,$(target)) : %.html : %.jemdoc pnpg.conf MENU
	jemdoc $(FLAGS) $<

