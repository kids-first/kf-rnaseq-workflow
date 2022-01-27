cwlVersion: v1.0
class: ExpressionTool
id: expression_get_rg_sm
requirements:
  - class: InlineJavascriptRequirement
inputs:
  rg:
    type: File
    inputBinding: {loadContents: true}
outputs:
  output: { type: 'string[]' }

expression: |
  ${
    var lines = inputs.rg.contents.split('\n');
    var rg_sms = [];
    for (var i = 0; i < lines.length; i++) {
      if (lines[i]) {
        var items = lines[i].split('\t');
        for (var j = 0; j < items.length; j++) {
          if (items[j].startsWith('SM')) {
            var sample = items[j].split(':')[1];
            rg_sms.push(sample);
          }
        }
      }
    }
    return {output: rg_sms};
  }
