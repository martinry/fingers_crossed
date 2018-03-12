from selenium import webdriver

def find_element(driver, xpath):
    for x in range(0, 20):  # try 4 times
        str_error = 'always an errorrrr'
        try:
            element = driver.find_element_by_xpath(xpath)
            str_error = ''
            print('Todo chido')
        except Exception as e:
            str_error = 'Errooor'
            print(str(e))
            pass

        if str_error != '':
            time.sleep(2)  # wait for 2 seconds before trying to fetch the data again
            print('Lets try again')
        else:
            return element
    print('Super slow connection... Sorry')
    
    
def retrieve_from_pfam(seq):
    driver = webdriver.Firefox()
    driver.get('https://www.ebi.ac.uk/Tools/pfa/pfamscan/')
    
    seq_box = find_element(driver, "//textarea[@id='sequence']")
    seq_box.sendKeys(seq)
    
    submit = find.element(driver, "//input[@name='submit']")
    submit.click()
    
    output_link = find.element(driver, "//a[@href[contains(.,'/out')]]")
    output_link.click()
    
    output = find.element(driver, "//pre")
    print(output.text)
    
retrieve_from_pfam('MAAVILPSTAAPSSLFPASQQKGHTQGGELVNELLTSWLRGLVTFEDVAVEFTQEEWALLDPAQRTLYRDVMLENCRNLASLGCRVNKPSLISQLEQDKKVVTEERGILPSTCPDLETLLKAKWLTPKKNVFRKEQSKGVKTERSHRGVKLNECNQCFKVFSTKSNLTQHKRIHTGEKPYDCSQCGKSFSSRSYLTIHKRIHNGEKPYECNHCGKAFSDPSSLRLHLRIHTGEKPYECNQCFHVFRTSCNLKSHKRIHTGENHHECNQCGKAFSTRSSLTGHNSIHTGEKPYECHDCGKTFRKSSYLTQHVRTHTGEKPYECNECGKSFSSSFSLTVHKRIHTGEKPYECSDCGKAFNNLSAVKKHLRTHTGEKPYECNHCGKSFTSNSYLSVHKRIHNRWI')