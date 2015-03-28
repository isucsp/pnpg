
enum MenuDescriptor{   
    MAC("Beam Hardening Correction", "MACPanel","MAC");

    private MenuDescriptor(String menu, String panelName, String abr){
        this.menu = menu;
        this.panelName = panelName;
        this.abr=abr;
    }
    
    /**
     * Method menu
     */
    public final String menu;
    public final String abr;
    
    /**
     * Description of method.
     */
    public final String panelName;
    
    @Override
    public String toString()
    {
        return menu;
    }

    static final int size =MenuDescriptor.values().length;
}
