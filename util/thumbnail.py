# -*- coding: utf-8 -*-
import logging

try:
    from indigo import Indigo, IndigoException
    from indigo_renderer import IndigoRenderer
    from indigo_inchi import IndigoInchi
    _indigo = Indigo()
    _renderer = IndigoRenderer(_indigo)
    _indigo_inchi = IndigoInchi(_indigo)

    def InChI2Thumbnail(inchi, output_format='png'):
        _indigo.setOption('render-output-format', output_format)
        _indigo.setOption('render-image-size', 250, 200)
        _indigo.setOption('render-margins', 10, 10)
        _indigo.setOption('render-coloring', True)
        _indigo.setOption('render-bond-length', 50.0)
        _indigo.setOption('render-stereo-style', 'none')
        _indigo.setOption('render-implicit-hydrogens-visible', False)
        _indigo.setOption('render-label-mode', 'terminal-hetero')

        try:
            indigo_mol = _indigo_inchi.loadMolecule(inchi)
            indigo_mol.aromatize()
            indigo_mol.layout()

            buf = _renderer.renderToBuffer(indigo_mol)
            return buf.tostring()
        except IndigoException as e:
            logging.warning("Cannot draw structure: %s" % str(e))
            return None

except ImportError:

    def InChI2Thumbnail(inchi, output_format='png'):
        logging.error('Indigo is not installed, cannot draw structures.')
        return None

if __name__ == '__main__':
    #inchi = 'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1'
    inchi = 'InChI=1S/C5H10N2O3/c6-3(5(9)10)1-2-4(7)8/h3H,1-2,6H2,(H2,7,8)(H,9,10)/t3-/m0/s1'
    thumb = InChI2Thumbnail(inchi, output_format='svg')
    if thumb is not None:
        fp = open('temp.html', 'w')
        fp.write(thumb)
        fp.close()
