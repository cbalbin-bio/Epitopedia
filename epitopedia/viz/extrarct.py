from gemmi import cif

doc = cif.read_file()
block = doc.sole_block()
new_doc = cif.Document()
new_block = new_doc.add_new_block(block.name)
new_block.add_item(block.find_loop_item("_atom_site.group_PDB"))

new_doc.as_string()
