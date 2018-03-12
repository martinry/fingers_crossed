from selenium import webdriver
import time, sys

def find_element(driver, xpath):
    for x in range(0, 3):  # try 4 times
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


def retrieve_from_ensembl(seq):
    driver = webdriver.Firefox()
    driver.get('http://www.ensembl.org/Multi/Search/Results?q=' + seq + ';site=ensembl_all')

    top_link = find_element(driver, "//a[@class='table_toplink']")
    top_link.click()

    orthologues = find_element(driver, "//a[@class='Orthologues']")
    orthologues.click()

    download_orthologues = find_element(driver, "//a[@class='export modal_link']")
    download_orthologues.click()

    fasta = find_element(driver, "//option[@value='FASTA']")
    fasta.click()

    preview = find_element(driver, "//input[@value='Preview']")
    preview.click()

    alignments = find_element(driver, "//pre[@style='color:#333']")
    with open('alignments_orthologs.fasta', 'w') as fout:
        fout.write(alignments.text)

    driver.close()
    
retrieve_from_ensembl("ENSP00000288602")
    