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
import { MatDialog } from '@angular/material/dialog';
import { Allotype } from '@app/shared/models/allotype/allotype.model';
import { AlleleExpansionComponent } from '../allele-expansion/allele-expansion.component';

@Component({
  selector: 'app-expression',
  templateUrl: './expression.component.html',
  styleUrls: ['./expression.component.scss']
})
export class ExpressionComponent implements OnInit {
  @Input() expression: string;
  @Input() type: string;
  @Input() transparent: boolean;
  @Input() movedPatientGenotype: boolean;
  @Input() allotype: Allotype;
  @Input() alleleDisplayed: string;

  constructor(public dialog: MatDialog) { }

  ngOnInit() {
  }

  openAllelePanel() {
    if (this.allotype.alleles){
      this.dialog.open(AlleleExpansionComponent, {
        data : {allotype : this.allotype},
        height: '80%'
      });
    }
  }

  processExpressionText(expressionText : string) : string{
    if (!expressionText){
      return null;
    }
    let str = expressionText;
    return str.replace('?', '').slice(0,4);
  }
}
