import qeschema

pw_document = qeschema.PwDocument()
pw_document.read('./aiida.xml')

dict_data = pw_document.to_dict()
print(dict_data['qes:espresso'])
