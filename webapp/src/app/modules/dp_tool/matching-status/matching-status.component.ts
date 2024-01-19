/*
 * Copyright (c) 2024 NMDP.
 *
 * This file is part of the NMDP DP Tool 
 * (see https://github.com/nmdp-bioinformatics/dpb1-expression).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
import { Component, OnInit, Input } from '@angular/core';
import { Subject } from '@app/shared/models/subject/subject.model';

@Component({
  selector: 'app-matching-status',
  templateUrl: './matching-status.component.html',
  styleUrls: ['./matching-status.component.scss']
})
export class MatchingStatusComponent implements OnInit {
  @Input() patient: Subject[];
  @Input() subject: Subject;
  tceIconPath : string;

  constructor() { }

  ngOnInit() {
  }

  getIconPath(mode : string){
    let filename = "";
    if (mode == 'expr'){
      const filenames = {'low' : 'expr-low',
                          'high' : 'expr-high',
                          'unknown' : 'unknown'}
      filename = filenames[this.subject.mismatchedAllelePatExprLevel];
    } else if (mode == 'tce'){
      const filenames = {'Permissive' : 'tce-perm',
                          'GvH_nonpermissive' : 'tce-NP',
                          'HvG_nonpermissive' : 'tce-NP'}
      filename = filenames[this.subject.tceMatch];
    }
  }

  getAllotypes(index : number){
    const donorAllotype = this.subject.allotypes.filter(a => a.sharedIndex != null)[0];
    const sharedIndex = donorAllotype.sharedIndex;
    if (index == 0){
      return [this.patient[0].allotypes[1 - sharedIndex]];
    } else if (index == 1){
      return this.subject.allotypes.filter(a => a.sharedIndex == null);
    } else {
      return this.subject.allotypes.filter(a => a.sharedIndex != null)
              .concat(this.patient[0].allotypes[sharedIndex]);
    }
  }

}
