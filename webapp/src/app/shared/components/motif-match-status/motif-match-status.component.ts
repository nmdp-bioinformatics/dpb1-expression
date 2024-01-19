import { Component, Input, OnInit } from '@angular/core';

@Component({
  selector: 'app-motif-match-status',
  templateUrl: './motif-match-status.component.html',
  styleUrls: ['./motif-match-status.component.scss']
})
export class MotifMatchStatusComponent implements OnInit {
  @Input() data: string;
  @Input() annotated: boolean;
  @Input() mode: string;
  basePath : string = "assets/img/";

  constructor() { }

  ngOnInit() {
  }

  getInfo() : Object{
    let filename = "";
    let filenames = {};
    let matchCategories = {};
    let matchCategoriesExpanded = {};
    if (this.mode == 'expr'){
      filenames = {'favorable' : 'expr-low',
                    'unfavorable' : 'expr-high'};
      matchCategories = {'favorable' : 'Expression favorable',
                        'unfavorable' : 'Expression unfavorable',
                        'unknown' : 'Expression unknown',
                        'mismatched' : 'Expression irrelevant (double mismatch)',
                        'matched' : 'Expression irrelevant (allele match)'};
    } else if (this.mode == 'tce'){
      filenames = {'Permissive' : 'tce-perm',
                    'GvH_nonpermissive' : 'tce-NP-GvH',
                    'HvG_nonpermissive' : 'tce-NP-HvG'};
      matchCategories = {'Permissive' : 'TCE permissive',
                    'GvH_nonpermissive' : 'TCE NP (GvH)',
                    'HvG_nonpermissive' : 'TCE NP (HvG)',
                    'Allele' : 'Matched'};
      matchCategoriesExpanded = {'Permissive' : 'permissive',
                    'GvH_nonpermissive' : 'non-permissive (graft-versus-host)',
                    'HvG_nonpermissive' : 'non-permissive (host-versus-graft)',
                    'Allele' : 'allele matched'};
    }
    filename = filenames[this.data];
    let filepath = null;
    if (filename){
      filepath = this.basePath + filename + '.png'
    }
    return {'filepath' : filepath,
     'category' : matchCategories[this.data],
     'categoryExpanded' : matchCategoriesExpanded[this.data]};
  }

}
